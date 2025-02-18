"""
This is a central script for the stages associated with Talos (previously AIP)

Due to privacy concerns, any configuration elements specific to an individual
project are stored in a separate config file, which is in a private repository.
 - production-pipelines-configuration/blob/main/config/talos/talos.toml

The cohort/project specific elements are described at the bottom of the Talos default config file here:
 - production-pipelines/blob/main/configs/default/talos.toml

[cohorts.DATASET_NAME]
- cohort_panels, list[int], PanelApp panel IDs to apply to this project. By default, only the Mendeliome (137) is
  applied to all analyses.

[cohorts.DATASET_NAME.genome]  # repeated for exome if appropriate
- DATASET_NAME is taken from config[workflow][dataset]
- historic_results, str, path to historic results directory for this project. If present, this folder will be used to
  detect the results from the latest run (so that variants which have been previously seen are annotated with the date
  they were first seen), and at the conclusion of the analysis a new version of the file is written with an additional
  new results from this run.
  If this config entry is absent, no previous results will be identified, so all results will have a 'first_seen'
  annotation of today's date, and the results from the current run won't be able to inform the next run.
- seqr_instance, str, URL for this seqr instance. Remove if not in seqr
- seqr_project, str, project ID for this project/seq type. Remove if not in seqr
- seqr_lookup, str, path to JSON file containing the mapping of CPG internal IDs to Seqr IDs, as generated by
  talos/helpers/process_seqr_metadata.py. If a new Seqr load has been performed, this file will need to be updated. This
  is stored in the test bucket so that it can be overwritten by a user. Remove if cohort is not in Seqr.

Takes as input:
    - Annotated MT path (either found in metamist or directly in config)
    - HPO.obo file (read from cpg-common reference location)
    - Seqr<->Internal ID mapping file (if appropriate, from Config)
    - History folder (if appropriate, from Config)
Generates:
    - Config file for this run
    - PED file for this cohort (extended format, with additional columns for Ext. ID and HPO terms)
    - Latest panels found through HPO matches against this cohort
    - PanelApp results
    - A mapping of all the Gene IDs in PanelApp to their Gene Symbol
    - Category-labelled VCF
    - Talos results JSON (metamist `aip-results` analysis)
    - Talos results + phenotype-match annotation
    - Talos report HTML (metamist `aip-report` analysis)

This will have to be run using Full permissions as we will need to reference data in test and main buckets.
"""

from datetime import datetime
from functools import lru_cache
from functools import cache
from os.path import join
from random import randint

import toml

from cloudpathlib.anypath import to_anypath

from cpg_utils import Path
from cpg_utils.config import ConfigError, config_retrieve
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch
from cpg_flow.resources import STANDARD
from cpg_flow.targets.dataset import Dataset
from cpg_flow.utils import exists, get_logger, tshirt_mt_sizing
from cpg_flow.stage import DatasetStage, StageInput, StageOutput, stage
from metamist.graphql import gql, query

CHUNKY_DATE = datetime.now().strftime('%Y-%m-%d')  # noqa: DTZ005
MTA_QUERY = gql(
    """
    query MyQuery($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
                meta
            }
        }
    }
""",
)
# used when building a runtime configuration
SEQR_KEYS: list[str] = ['seqr_project', 'seqr_instance', 'seqr_lookup']


@lru_cache(maxsize=1)
def get_date_string() -> str:
    """
    allows override of the date folder to continue/re-run previous analyses

    Returns:
        either an override in config, or the default (today, YYYY-MM-DD)
    """
    return config_retrieve(['workflow', 'date_folder_override'], CHUNKY_DATE)


@lru_cache(1)
def get_date_folder() -> str:
    """
    allows override of the date folder to continue/re-run previous analyses
    Returns:
        either an override in config, or the default "reanalysis/(today, YYYY-MM-DD)"
    """
    return join('reanalysis', get_date_string())


@cache
def query_for_sv_mt(dataset: str) -> tuple[str, str] | None:
    """
    query for the latest SV MT for a dataset
    bonus - is we're searching for CNVs, we search for multiple
    return the full paths and filenames only, as 2 lists

    Args:
        dataset (str): project to query for

    Returns:
        the path & filename to the latest MT for the given type, or None
    """

    # this separation is based on the gCNV(exome)/GATK-SV(genome) divide at CPG
    sv_type = 'cnv' if config_retrieve(['workflow', 'sequencing_type']) == 'exome' else 'sv'

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    # we only want the final stage MT, subset to the specific dataset
    final_stage_lookup = {'cnv': 'AnnotateDatasetCNV', 'sv': 'AnnotateDatasetSv'}

    result = query(MTA_QUERY, variables={'dataset': query_dataset, 'type': sv_type})
    mt_by_date: dict[str, str] = {}
    for analysis in result['project']['analyses']:
        if (
            analysis['output']
            and analysis['output'].endswith('.mt')
            and (analysis['meta']['sequencing_type'] == config_retrieve(['workflow', 'sequencing_type']))
            and (analysis['meta']['stage'] == final_stage_lookup[sv_type])
        ):
            mt_by_date[analysis['timestampCompleted']] = analysis['output']

    # perfectly acceptable to not have an input SV MT
    if not mt_by_date:
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    sv_file = mt_by_date[sorted(mt_by_date)[-1]]
    filename = sv_file.split('/')[-1]
    return sv_file, filename


@cache
def query_for_latest_hail_object(
    dataset: str,
    analysis_type: str,
    object_suffix: str = '.mt',
    exact_string: str | None = None,
) -> str:
    """
    query for the latest MT for a dataset
    the exact metamist entry type to search for is handled by config, defaulting to the new rd_combiner MT
    Args:
        dataset (str):       project to query for
        analysis_type (str): analysis type to query for - rd_combiner writes MTs to metamist as 'matrixtable',
                             seqr_loader used 'custom': using a config entry we can decide which type to use
        object_suffix (str): suffix to detect on analysis path (choose between .mt or .ht)
        exact_string (str):  if set, analysis entries will only be considered if this substring is in the path
                             the intended use here is targeting outputs from a single Stage (annotate cohort vs dataset)
    Returns:
        str, the path to the latest object for the given type
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset
    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    get_logger().info(f'Querying for {analysis_type} in {query_dataset}')

    result = query(MTA_QUERY, variables={'dataset': query_dataset, 'type': analysis_type})

    # get all the relevant entries, and bin by date
    mt_by_date = {}
    for analysis in result['project']['analyses']:
        # if exact string is absent, or set and in this output path, use it
        if (
            analysis['output']
            and analysis['output'].endswith(object_suffix)
            and (analysis['meta']['sequencing_type'] == config_retrieve(['workflow', 'sequencing_type']))
            and (exact_string is None or exact_string in analysis['output'])
        ):
            mt_by_date[analysis['timestampCompleted']] = analysis['output']

    if not mt_by_date:
        raise ValueError(f'No MT found for dataset {query_dataset}')

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return mt_by_date[sorted(mt_by_date)[-1]]


@lru_cache(2)
def get_clinvar_table(key: str = 'clinvar_decisions') -> str:
    """
    this is used to retrieve two types of object - clinvar_decisions & clinvar_pm5

    try and identify the clinvar table to use
    - try the config specified path
    - fall back to storage:common default path - try with multiple dir names in case we change this
    - if neither works, choose to fail instead

    Args
        key (str): the key to look for in the config

    Returns:
        a path to a clinvar table, or None
    """

    if (clinvar_table := config_retrieve(['workflow', key], None)) is not None:
        get_logger().info(f'Using clinvar table {clinvar_table} from config')
        return clinvar_table

    get_logger().info(f'No forced {key} table available, trying default')

    # try multiple path variations - legacy dir name is 'aip_clinvar', but this may also change
    for default_name in ['clinvarbitration', 'aip_clinvar']:
        clinvar_table = join(
            config_retrieve(['storage', 'common', 'analysis']),
            default_name,
            datetime.now().strftime('%y-%m'),  # noqa: DTZ005
            f'{key}.ht',
        )

        if to_anypath(clinvar_table).exists():
            get_logger().info(f'Using clinvar table {clinvar_table}')
            return clinvar_table

    raise ValueError('no Clinvar Tables were identified')


@stage
class GeneratePED(DatasetStage):
    """
    revert to just using the metamist/CPG-flow Pedigree generation
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'pedigree.ped'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        expected_out = self.expected_outputs(dataset)
        pedigree = dataset.write_ped_file(out_path=expected_out)
        get_logger().info(f'PED file for {dataset.name} written to {pedigree}')

        return self.make_outputs(dataset, data=expected_out)


@stage
class MakeRuntimeConfig(DatasetStage):
    """
    create a config file for this run,
    this new config should include all elements specific to this Project and sequencing_type
    this new unambiguous config file should be used in all jobs
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'config.toml'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        # start off with a fresh config dictionary, including generic content
        new_config: dict = {
            'categories': config_retrieve(['categories']),
            'GeneratePanelData': config_retrieve(['GeneratePanelData']),
            'FindGeneSymbolMap': config_retrieve(['FindGeneSymbolMap']),
            'RunHailFiltering': config_retrieve(['RunHailFiltering']),
            'ValidateMOI': config_retrieve(['ValidateMOI']),
            'HPOFlagging': config_retrieve(['HPOFlagging']),
            'CreateTalosHTML': {},
        }

        # pull the content relevant to this cohort + sequencing type (mandatory in CPG)
        seq_type = config_retrieve(['workflow', 'sequencing_type'])
        dataset_conf = config_retrieve(['cohorts', dataset.name])
        seq_type_conf = dataset_conf.get(seq_type, {})

        # forbidden genes and forced panels
        new_config['GeneratePanelData']['forbidden_genes'] = dataset_conf.get('forbidden', [])
        new_config['GeneratePanelData']['forced_panels'] = dataset_conf.get('forced_panels', [])
        new_config['GeneratePanelData']['blacklist'] = dataset_conf.get('blacklist', None)

        # optionally, all SG IDs to remove from analysis
        new_config['ValidateMOI']['solved_cases'] = dataset_conf.get('solved_cases', [])

        # these attributes are present, or missing completely
        if all(x in seq_type_conf for x in SEQR_KEYS):
            for key in SEQR_KEYS:
                if key in seq_type_conf:
                    new_config['CreateTalosHTML'][key] = seq_type_conf[key]
        if 'external_labels' in seq_type_conf:
            new_config['CreateTalosHTML']['external_labels'] = seq_type_conf['external_labels']

        # add a location for the run history files
        if 'result_history' in seq_type_conf:
            new_config['result_history'] = seq_type_conf['result_history']

        expected_outputs = self.expected_outputs(dataset)

        with expected_outputs.open('w') as write_handle:
            toml.dump(new_config, write_handle)

        return self.make_outputs(target=dataset, data=expected_outputs)


@stage
class MakePhenopackets(DatasetStage):
    """
    this calls the script which reads phenotype data from metamist
    and generates a phenopacket file (GA4GH compliant)
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'phenopackets.json'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        generate a pedigree from metamist
        script to generate an extended pedigree format - additional columns for Ext. ID and HPO terms
        """
        job = get_batch().new_job('Generate Phenopackets from Metamist')
        job.cpu(1).image(config_retrieve(['workflow', 'driver_image']))

        expected_out = self.expected_outputs(dataset)
        query_dataset = dataset.name
        if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
            query_dataset += '-test'

        hpo_file = get_batch().read_input(config_retrieve(['GeneratePanelData', 'obo_file']))

        # mandatory argument
        seq_type = config_retrieve(['workflow', 'sequencing_type'])

        # insert a little stagger
        job.command(f'sleep {randint(0, 30)}')

        job.command(
            f'MakePhenopackets --dataset {query_dataset} --output {job.output} --type {seq_type} --hpo {hpo_file}',
        )
        get_batch().write_output(job.output, str(expected_out))
        get_logger().info(f'Phenopacket file for {dataset.name} going to {expected_out}')

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[MakePhenopackets, MakeRuntimeConfig])
class GeneratePanelData(DatasetStage):
    """
    PythonJob to find HPO-matched panels
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        """
        only one output, the panel data
        """
        return dataset.prefix() / get_date_folder() / 'hpo_panel_data.json'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Find HPO-matched Panels: {dataset.name}')
        job.cpu(1).image(config_retrieve(['workflow', 'driver_image']))

        # use the new config file
        conf_in_batch = get_batch().read_input(inputs.as_str(dataset, MakeRuntimeConfig))

        expected_out = self.expected_outputs(dataset)

        hpo_file = get_batch().read_input(config_retrieve(['GeneratePanelData', 'obo_file']))
        local_phenopacket = get_batch().read_input(inputs.as_str(target=dataset, stage=MakePhenopackets))

        job.command(f'export TALOS_CONFIG={conf_in_batch}')
        # insert a little stagger

        job.command(f'sleep {randint(0, 30)}')
        job.command(f'GeneratePanelData --input {local_phenopacket} --output {job.output} --hpo {hpo_file}')
        get_batch().write_output(job.output, str(expected_out))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[GeneratePanelData, MakeRuntimeConfig])
class QueryPanelapp(DatasetStage):
    """
    query PanelApp for up-to-date gene lists
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'panelapp_data.json'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Query PanelApp: {dataset.name}')
        job.cpu(1).image(config_retrieve(['workflow', 'driver_image']))

        # use the new config file
        runtime_config = inputs.as_str(dataset, MakeRuntimeConfig)
        conf_in_batch = get_batch().read_input(runtime_config)

        # read the previous Stage's output into the batch
        hpo_panel_json = get_batch().read_input(inputs.as_str(target=dataset, stage=GeneratePanelData))
        expected_out = self.expected_outputs(dataset)
        job.command(f'export TALOS_CONFIG={conf_in_batch}')

        # insert a little stagger
        job.command(f'sleep {randint(20, 300)}')
        job.command(f'QueryPanelapp --input {hpo_panel_json} --output {job.output}')
        get_batch().write_output(job.output, str(expected_out))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[MakeRuntimeConfig, QueryPanelapp])
class FindGeneSymbolMap(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'symbol_to_ensg.json'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        expected_out = self.expected_outputs(dataset)

        job = get_batch().new_job(f'Find Symbol-ENSG lookup: {dataset.name}')
        job.cpu(1).image(config_retrieve(['workflow', 'driver_image']))

        # use the new config file
        conf_in_batch = get_batch().read_input(inputs.as_str(dataset, MakeRuntimeConfig))

        panel_json = inputs.as_str(target=dataset, stage=QueryPanelapp)
        job.command(f'export TALOS_CONFIG={conf_in_batch}')

        # insert a little stagger
        job.command(f'sleep {randint(0, 30)}')
        job.command(f'FindGeneSymbolMap --input {panel_json} --output {job.output}')

        get_batch().write_output(job.output, str(expected_out))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[QueryPanelapp, GeneratePED, MakeRuntimeConfig])
class RunHailFiltering(DatasetStage):
    """
    hail job to filter & label the MT
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'hail_labelled.vcf.bgz'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        input_mt = config_retrieve(
            ['workflow', 'matrix_table'],
            default=query_for_latest_hail_object(
                dataset=dataset.name,
                analysis_type=config_retrieve(['workflow', 'mt_entry_type'], default='matrixtable'),
                object_suffix='.mt',
            ),
        )

        job = get_batch().new_job(f'Run hail labelling: {dataset.name}')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command('set -eux pipefail')

        # time in seconds before this jobs self-destructs
        # some recent runs of this have hit the GCP copy and... hung indefinitely at great expense
        # the highest current runtime when successful is just shy of 4 hours
        job.timeout(config_retrieve(['RunHailFiltering', 'timeouts', 'small_variants'], 15000))

        # use the new config file
        runtime_config = inputs.as_str(dataset, MakeRuntimeConfig)
        conf_in_batch = get_batch().read_input(runtime_config)

        # MTs can vary from <10GB for a small exome, to 170GB for a larger one, Genomes ~500GB
        required_storage = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(dataset.get_sequencing_group_ids()),
        )
        required_cpu: int = config_retrieve(['RunHailFiltering', 'cores', 'small_variants'], 8)
        required_mem: str = config_retrieve(['RunHailFiltering', 'memory', 'small_variants'], 'highmem')

        # doubling storage due to the repartitioning
        job.cpu(required_cpu).storage(f'{required_storage * 2}Gi').memory(required_mem)

        panelapp_json = get_batch().read_input(inputs.as_str(target=dataset, stage=QueryPanelapp))
        pedigree = get_batch().read_input(inputs.as_str(target=dataset, stage=GeneratePED))
        expected_out = self.expected_outputs(dataset)

        # copy vcf & index out manually
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # find the clinvar tables, and localise
        clinvar_decisions = get_clinvar_table()
        clinvar_name = clinvar_decisions.split('/')[-1]

        # localise the clinvar decisions table
        job.command(f'gcloud --no-user-output-enabled storage cp -r {clinvar_decisions} $BATCH_TMPDIR')
        job.command('echo "ClinvArbitration decisions copied"')

        # see if we can find any exomiser results to integrate
        try:
            exomiser_ht = query_for_latest_hail_object(
                dataset=dataset.name,
                analysis_type='exomiser',
                object_suffix='.ht',
            )
            exomiser_name = exomiser_ht.split('/')[-1]
            job.command(f'gcloud --no-user-output-enabled storage cp -r {exomiser_ht} $BATCH_TMPDIR')
            job.command('echo "Exomiser HT copied"')
            exomiser_argument = f'--exomiser "${{BATCH_TMPDIR}}/{exomiser_name}" '
        except ValueError:
            get_logger().info(f'No exomiser results found for {dataset.name}, skipping')
            exomiser_argument = ' '

        # find, localise, and use the SpliceVarDB table, if available - if not, don't pass the flag
        # currently just passed in from config, will eventually be generated a different way
        if svdb := config_retrieve(['RunHailFiltering', 'svdb_ht'], None):
            if not exists(svdb):
                raise ValueError(f'SVDB {svdb} does not exist')

            svdb_name = svdb.split('/')[-1]
            job.command(f'gcloud --no-user-output-enabled storage cp -r {svdb} $BATCH_TMPDIR')
            job.command('echo "SpliceVarDB MT copied"')
            svdb_argument = f'--svdb "${{BATCH_TMPDIR}}/{svdb_name}" '
        else:
            svdb_argument = ' '

        pm5 = get_clinvar_table('clinvar_pm5')
        pm5_name = pm5.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {pm5} $BATCH_TMPDIR')
        job.command('echo "ClinvArbitration PM5 copied"')

        # finally, localise the whole MT (this takes the longest
        mt_name = input_mt.split('/')[-1]
        job.command(f'gcloud --no-user-output-enabled storage cp -r {input_mt} $BATCH_TMPDIR')
        job.command('echo "Cohort MT copied"')

        job.command(f'export TALOS_CONFIG={conf_in_batch}')
        job.command(
            'RunHailFiltering '
            f'--input "${{BATCH_TMPDIR}}/{mt_name}" '
            f'--panelapp {panelapp_json} '
            f'--pedigree {pedigree} '
            f'--output {job.output["vcf.bgz"]} '
            f'--clinvar "${{BATCH_TMPDIR}}/{clinvar_name}" '
            f'--pm5 "${{BATCH_TMPDIR}}/{pm5_name}" '
            f'--checkpoint "${{BATCH_TMPDIR}}/checkpoint.mt" '
            f'{svdb_argument} '
            f'{exomiser_argument} ',
        )
        get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[QueryPanelapp, GeneratePED, MakeRuntimeConfig])
class RunHailFilteringSV(DatasetStage):
    """
    hail job to filter & label the SV MT
    """

    def expected_outputs(self, dataset: Dataset) -> Path | None:
        sv_input = query_for_sv_mt(dataset.name)
        if sv_input:
            _filepath, filename = sv_input
            return dataset.prefix() / get_date_folder() / f'label_{filename}.vcf.bgz'
        return None

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        expected_out = self.expected_outputs(dataset)

        # nothing to do, no input SV data
        if not expected_out:
            return self.make_outputs(dataset, skipped=True)

        conf_in_batch = get_batch().read_input(inputs.as_str(dataset, MakeRuntimeConfig))
        panelapp_json = get_batch().read_input(inputs.as_str(target=dataset, stage=QueryPanelapp))
        local_ped = get_batch().read_input(inputs.as_str(target=dataset, stage=GeneratePED))

        required_storage: int = config_retrieve(['RunHailFiltering', 'storage', 'sv'], 10)
        required_cpu: int = config_retrieve(['RunHailFiltering', 'cores', 'sv'], 2)

        # query for (cached) SV results
        sv_results = query_for_sv_mt(dataset.name)

        # if there were no SV results, expected_out was empty so we returned above
        assert sv_results, f'No SV results found for {dataset.name}'
        sv_path, sv_file = sv_results
        job = get_batch().new_job(f'Run hail SV labelling: {dataset.name}, {sv_file}')
        # manually extract the VCF and index
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        job.image(config_retrieve(['workflow', 'driver_image']))
        # generally runtime under 10 minutes
        job.timeout(config_retrieve(['RunHailFiltering', 'timeouts', 'sv'], 3600))

        # use the new config file
        STANDARD.set_resources(job, ncpu=required_cpu, storage_gb=required_storage, mem_gb=16)

        # copy the mt in
        job.command(f'gcloud --no-user-output-enabled storage cp -r {sv_path} $BATCH_TMPDIR')
        job.command(f'export TALOS_CONFIG={conf_in_batch}')
        job.command(
            'RunHailFilteringSV '
            f'--input "${{BATCH_TMPDIR}}/{sv_file}" '
            f'--panelapp {panelapp_json} '
            f'--pedigree {local_ped} '
            f'--output {job.output["vcf.bgz"]} ',
        )
        get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[
        GeneratePED,
        GeneratePanelData,
        QueryPanelapp,
        RunHailFiltering,
        RunHailFilteringSV,
        MakeRuntimeConfig,
    ],
)
class ValidateMOI(DatasetStage):
    """
    run the labelled VCF -> results JSON stage
    """

    def expected_outputs(self, dataset: Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'summary_output.json'

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Talos summary: {dataset.name}')
        job.cpu(config_retrieve(['talos_stages', 'ValidateMOI', 'cpu'], 2.0)).memory(
            config_retrieve(['talos_stages', 'ValidateMOI', 'memory'], 'highmem'),
        ).storage(config_retrieve(['talos_stages', 'ValidateMOI', 'storage'], '10Gi')).image(
            config_retrieve(['workflow', 'driver_image']),
        )
        # use the new config file
        runtime_config = inputs.as_str(dataset, MakeRuntimeConfig)
        conf_in_batch = get_batch().read_input(runtime_config)

        hpo_panels = get_batch().read_input(inputs.as_str(dataset, GeneratePanelData))
        pedigree = get_batch().read_input(inputs.as_str(target=dataset, stage=GeneratePED))
        hail_inputs = inputs.as_str(dataset, RunHailFiltering)

        # If there are SV VCFs, read each one in and add to the arguments
        if query_for_sv_mt(dataset.name) is None:
            get_logger().warning(f'No SV results found for {dataset.name}')
            sv_vcf_arg = ''
        else:
            # only go looking for inputs from prior stage where we expect to find them
            hail_sv_input = inputs.as_str(dataset, RunHailFilteringSV)
            labelled_sv_vcf = get_batch().read_input_group(
                **{'vcf.bgz': hail_sv_input, 'vcf.bgz.tbi': f'{hail_sv_input}.tbi'},
            )['vcf.bgz']

            sv_vcf_arg = f'--labelled_sv {labelled_sv_vcf}'

        panel_input = get_batch().read_input(inputs.as_str(dataset, QueryPanelapp))
        labelled_vcf = get_batch().read_input_group(
            **{
                'vcf.bgz': hail_inputs,
                'vcf.bgz.tbi': hail_inputs + '.tbi',
            },
        )['vcf.bgz']

        job.command(f'export TALOS_CONFIG={conf_in_batch}')
        job.command(
            'ValidateMOI '
            f'--labelled_vcf {labelled_vcf} '
            f'--output {job.output} '
            f'--panelapp {panel_input} '
            f'--pedigree {pedigree} '
            f'--participant_panels {hpo_panels} '
            f'{sv_vcf_arg}',
        )
        expected_out = self.expected_outputs(dataset)
        get_batch().write_output(job.output, str(expected_out))
        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[MakeRuntimeConfig, FindGeneSymbolMap, ValidateMOI],
    analysis_type='aip-results',
    analysis_keys=['pheno_annotated', 'pheno_filtered'],
)
class HPOFlagging(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        date_prefix = dataset.prefix() / get_date_folder()
        return {
            'pheno_annotated': date_prefix / 'pheno_annotated_report.json',
            'pheno_filtered': date_prefix / 'pheno_filtered_report.json',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(dataset)

        # TODO IDK, get these from config?
        phenio_db = get_batch().read_input(config_retrieve(['HPOFlagging', 'phenio_db']))
        gene_to_phenotype = get_batch().read_input(config_retrieve(['HPOFlagging', 'gene_to_phenotype']))

        job = get_batch().new_job(f'Label phenotype matches: {dataset.name}')
        job.cpu(2.0).memory('highmem').image(config_retrieve(['workflow', 'driver_image'])).storage('20Gi')

        # use the new config file
        runtime_config = inputs.as_str(dataset, MakeRuntimeConfig)
        conf_in_batch = get_batch().read_input(runtime_config)

        results_json = get_batch().read_input(inputs.as_str(dataset, ValidateMOI))
        gene_map = get_batch().read_input(inputs.as_str(dataset, FindGeneSymbolMap))

        job.command(f'export TALOS_CONFIG={conf_in_batch}')
        job.command(
            'HPOFlagging '
            f'--input {results_json} '
            f'--gene_map {gene_map} '
            f'--gen2phen {gene_to_phenotype} '
            f'--phenio {phenio_db} '
            f'--output {job.output} '
            f'--phenout {job.phenout} ',
        )

        get_batch().write_output(job.output, str(outputs['pheno_annotated']))
        get_batch().write_output(job.phenout, str(outputs['pheno_filtered']))

        return self.make_outputs(target=dataset, jobs=job, data=outputs)


@stage(
    required_stages=[HPOFlagging, QueryPanelapp, MakeRuntimeConfig],
    analysis_type='aip-report',
    analysis_keys=['results_html', 'latest_html'],
    tolerate_missing_output=True,
)
class CreateTalosHTML(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        date_folder_prefix = dataset.prefix(category='web') / get_date_folder()
        return {
            'results_html': date_folder_prefix / 'summary_output.html',
            'latest_html': date_folder_prefix / f'summary_latest_{get_date_string()}.html',
            'folder': date_folder_prefix,
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Talos HTML: {dataset.name}')
        job.image(config_retrieve(['workflow', 'driver_image'])).memory('standard').cpu(1.0)

        # use the new config file
        runtime_config = inputs.as_str(dataset, MakeRuntimeConfig)
        conf_in_batch = get_batch().read_input(runtime_config)

        results_json = get_batch().read_input(str(inputs.as_dict(dataset, HPOFlagging)['pheno_annotated']))
        panel_input = get_batch().read_input(inputs.as_str(dataset, QueryPanelapp))
        expected_out = self.expected_outputs(dataset)

        # this + copy_common_env (called by default) will be enough to do a gcloud copy of the outputs
        authenticate_cloud_credentials_in_job(job)

        # this will write output files directly to GCP
        # report splitting is arbitrary, so can't be reliably done in Hail
        job.command(f'export TALOS_CONFIG={conf_in_batch}')

        # create a new directory for the results
        job.command('mkdir html_outputs')
        job.command('cd html_outputs')

        command_string = (
            'CreateTalosHTML '
            f'--input {results_json} '
            f'--panelapp {panel_input} '
            f'--output {expected_out["results_html"].name} '
            f'--latest {expected_out["latest_html"].name} '
        )

        if report_splitting := config_retrieve(['workflow', 'report_splitting', dataset.name], False):
            command_string += f' --split_samples {report_splitting}'

        job.command(command_string)

        # copy the results to the bucket
        expected_out_folder_string = str(expected_out['folder'])
        job.command(f'gcloud storage cp -r * {expected_out_folder_string}')

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[ValidateMOI, MakeRuntimeConfig],
    analysis_keys=['seqr_file', 'seqr_pheno_file'],
    analysis_type='custom',
    tolerate_missing_output=True,
)
class MinimiseOutputForSeqr(DatasetStage):
    """
    takes the results file from Seqr and produces a minimised form for Seqr ingestion
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        analysis_folder_prefix = dataset.prefix(category='analysis') / 'seqr_files'
        return {
            'seqr_file': analysis_folder_prefix / f'{get_date_folder()}_seqr.json',
            'seqr_pheno_file': analysis_folder_prefix / f'{get_date_folder()}_seqr_pheno.json',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        # pull out the config section relevant to this datatype & cohort
        # if it doesn't exist for this sequencing type, fail gracefully
        seq_type = config_retrieve(['workflow', 'sequencing_type'])
        try:
            seqr_lookup = config_retrieve(['cohorts', dataset.name, seq_type, 'seqr_lookup'])
        except ConfigError:
            get_logger().warning(f'No Seqr lookup file for {dataset.name} {seq_type}')
            return self.make_outputs(dataset, skipped=True)

        input_localised = get_batch().read_input(inputs.as_str(dataset, ValidateMOI))

        # create a job to run the minimisation script
        job = get_batch().new_job(f'Talos Prep for Seqr: {dataset.name}')
        job.image(config_retrieve(['workflow', 'driver_image'])).cpu(1.0).memory('lowmem')

        # use the new config file
        runtime_config = inputs.as_str(dataset, MakeRuntimeConfig)
        conf_in_batch = get_batch().read_input(runtime_config)

        lookup_in_batch = get_batch().read_input(seqr_lookup)
        job.command(f'export TALOS_CONFIG={conf_in_batch}')
        job.command(
            'MinimiseOutputForSeqr '
            f'--input {input_localised} '
            f'--output {job.out_json} '
            f'--pheno {job.pheno_json} '
            f'--external_map {lookup_in_batch}',
        )

        # write the results out
        expected_out = self.expected_outputs(dataset)
        get_batch().write_output(job.out_json, str(expected_out['seqr_file']))
        get_batch().write_output(job.pheno_json, str(expected_out['seqr_pheno_file']))
        return self.make_outputs(dataset, data=expected_out, jobs=job)
