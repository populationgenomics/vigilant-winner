"""
Hail Batch jobs to run VEP on a VCF in parallel.
"""

import logging
from enum import Enum
from typing import Literal, Optional, Union, Dict, Tuple, List

from cloudpathlib import CloudPath

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    image_path,
    reference_path,
    authenticate_cloud_credentials_in_job,
    query_command,
    fasta_res_group,
)

from . import query

logger = logging.getLogger(__file__)


class SequencingType(Enum):
    """
    Type (scope) of a sequencing experiment.
    """

    GENOME = 'genome'
    EXOME = 'exome'


def vep_jobs(  # pylint: disable=too-many-arguments
    b: Batch,
    vcf_path: CloudPath,
    hail_billing_project: str,
    hail_bucket: CloudPath,
    tmp_bucket: CloudPath,
    out_path: Optional[CloudPath] = None,
    overwrite: bool = False,
    scatter_count: Optional[int] = 50,
    sequencing_type: SequencingType = SequencingType.GENOME,
    intervals_path: Optional[CloudPath] = None,
    job_attrs: Optional[dict] = None,
) -> List[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """
    to_hail_table = out_path and out_path.suffix == '.ht'
    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    if out_path and not overwrite and out_path.exists():
        return [b.new_job('VEP [reuse]', job_attrs)]

    jobs: List[Job] = []
    intervals_j, intervals = get_intervals(
        b=b,
        cache_bucket=reference_path('intervals_prefix')
        / sequencing_type.value
        / f'{scatter_count}intervals',
        sequencing_type=sequencing_type,
        intervals_path=intervals_path,
        scatter_count=scatter_count,
    )
    jobs.append(intervals_j)

    vcf = b.read_input_group(
        **{'vcf.gz': str(vcf_path), 'vcf.gz.tbi': str(vcf_path) + '.tbi'}
    )

    parts_bucket = tmp_bucket / 'vep' / 'parts'
    part_files = []

    # Splitting variant calling by intervals
    for idx in range(scatter_count):

        # find the eventual output path if appropriate
        if to_hail_table:
            part_path = parts_bucket / f'part{idx + 1}.json_list'
        else:
            part_path = None

        # here we assume that if the eventual path exists, the subset and
        # annotation were both done and can be re-used
        # the subset-vcf is not persisted, so we can skip either both jobs,
        # or neither
        if part_path and part_path.exists() and not overwrite:
            part_files.append(part_path)
            continue

        subset_j = subset_vcf(
            b,
            vcf=vcf,
            interval=intervals[idx],
            job_attrs=job_attrs or dict(part=f'{idx + 1}/{scatter_count}'),
        )
        jobs.append(subset_j)
        # noinspection PyTypeChecker
        j = vep_one(
            b,
            vcf=subset_j.output_vcf['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=part_path,
            job_attrs=job_attrs or dict(part=f'{idx + 1}/{scatter_count}'),
            overwrite=overwrite,
        )
        jobs.append(j)
        if to_hail_table:
            part_files.append(part_path)
        else:
            part_files.append(j.output['vcf.gz'])

    if to_hail_table:

        # if already generated, don't regen
        if out_path and CloudPath(out_path).exists() and not overwrite:
            return jobs

        gather_j = gather_vep_json_to_ht(
            b=b,
            vep_results_paths=part_files,
            hail_billing_project=hail_billing_project,
            hail_bucket=hail_bucket,
            out_path=out_path,
            job_attrs=job_attrs,
        )
    else:
        assert len(part_files) == scatter_count
        gather_j, _gather_vcf = gather_vcfs(
            b=b,
            input_vcfs=part_files,
            out_vcf_path=out_path,
        )
    gather_j.depends_on(*jobs)
    jobs.append(gather_j)
    return jobs


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    intervals_path: Optional[CloudPath] = None,
    cache_bucket: Optional[CloudPath] = None,
    sequencing_type: SequencingType = SequencingType.GENOME,
    job_attrs: Optional[dict] = None,
) -> Tuple[Job, List[hb.Resource]]:
    """
    Add a job that split genome into partitions for variant calling parallelisation.

    Takes `intervals_path` if provided, otherwise calls `reference_path()`
    for the intervals of provided `sequencing_type`.

    Caches intervals for each partition into `cache_bucket`, if provided.

    This job calls picard's IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob
    /bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredicted number of intervals. WDL can handle
    that, but Hail Batch is not dynamic and have to expect certain number of output
    files.
    """
    job_attrs = job_attrs or dict(tool='picard_IntervalListTools')
    j = b.new_job(f'Make {scatter_count} intervals', job_attrs)

    if cache_bucket:
        # Checking previously cached split intervals.
        if (cache_bucket / '1.interval_list').exists():
            j.name += ' [use cached]'
            return j, [
                b.read_input(str(cache_bucket / f'{idx + 1}.interval_list'))
                for idx in range(scatter_count)
            ]

    # Taking intervals file for the sequencing_type.
    intervals_path = intervals_path or reference_path(
        f'broad/{sequencing_type.value}_calling_interval_lists',
    )

    j.image(image_path('picard'))
    j.memory('16Gi')
    j.storage('50G')
    j.cpu(4)

    break_bands_at_multiples_of = {
        SequencingType.GENOME: 100000,
        SequencingType.EXOME: 0,
    }.get(sequencing_type, 0)

    cmd = f"""
    mkdir /io/batch/out

    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(str(intervals_path))} \
    OUTPUT=/io/batch/out
    ls /io/batch/out
    ls /io/batch/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln /io/batch/out/{name}/scattered.interval_list {j[f'{idx + 1}.interval_list']}
        """

    j.command(cmd)
    if cache_bucket:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx + 1}.interval_list'],
                str(cache_bucket / f'{idx + 1}.interval_list'),
            )
    return j, [j[f'{idx + 1}.interval_list'] for idx in range(scatter_count)]


def subset_vcf(
    b: hb.Batch,
    vcf: hb.ResourceGroup,
    interval: hb.Resource,
    job_attrs: Optional[Dict] = None,
    output_vcf_path: Optional[CloudPath] = None,
) -> Job:
    """
    Subset VCF to provided intervals, and drop sample/genotype
    information (e.g. outputs sites-only VCF).
    """
    job_name = 'VEP: subset VCF'
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    j.memory('16Gi')
    j.storage('50G')
    j.cpu(2)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
    gatk SelectVariants \\
    -R {fasta_res_group(b).base} \\
    -V {vcf['vcf.gz']} \\
    -L {interval} \\
    -O /io/batch/tmp.vcf.gz

    gatk MakeSitesOnlyVcf \\
    -I /io/batch/tmp.vcf.gz \\
    -O {j.output_vcf['vcf.gz']} \\
    --CREATE_INDEX
    """
    j.command(cmd)
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j


def gather_vcfs(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceFile],
    overwrite: bool = True,
    out_vcf_path: Optional[CloudPath] = None,
    site_only: bool = False,
    gvcf_count: Optional[int] = None,
    job_attrs: Optional[dict] = None,
) -> Tuple[Job, hb.ResourceGroup]:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    if out_vcf_path and out_vcf_path.exists() and not overwrite:
        j.name += ' [reuse]'
        return j, b.read_input_group(
            **{
                'vcf.gz': str(out_vcf_path),
                'vcf.gz.tbi': f'{out_vcf_path}.tbi',
            }
        )

    j.memory('16Gi')
    j.storage('50G')
    j.cpu(16)
    if gvcf_count:
        j.storage((1 if site_only else 2) * gvcf_count)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v}' for v in input_vcfs])
    cmd = f"""
    # --ignore-safety-checks makes a big performance difference so we include it in
    # our invocation. This argument disables expensive checks that the file headers
    # contain the same set of genotyped samples and that files are in order
    # by position of first record.
    gatk --java-options -Xms25g \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
    {input_cmdl} \\
    --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    """
    j.command(cmd)
    if out_vcf_path:
        b.write_output(j.output_vcf, str(out_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def gather_vep_json_to_ht(
    b: Batch,
    vep_results_paths: List[CloudPath],
    hail_billing_project: str,
    hail_bucket: CloudPath,
    out_path: CloudPath,
    job_attrs: Optional[dict] = None,
):
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    j = b.new_job('VEP json to Hail table', job_attrs)
    j.image(get_config()['workflow']['driver_image'])
    cmd = query_command(
        query,
        query.vep_json_to_ht.__name__,
        [str(p) for p in vep_results_paths],
        str(out_path),
        setup_gcp=True,
        hail_billing_project=hail_billing_project,
        hail_bucket=str(hail_bucket),
    )
    j.command(cmd)
    return j


def vep_one(
    b: Batch,
    vcf: Union[CloudPath, hb.Resource],
    out_path: Optional[CloudPath] = None,
    out_format: Literal['vcf', 'json'] = 'vcf',
    job_attrs: Optional[Dict] = None,
    overwrite: bool = False,
) -> Job:
    """
    Run a single VEP job.
    """
    j = b.new_job('VEP', job_attrs)
    if out_path and out_path.exists() and not overwrite:
        j.name += ' [reuse]'
        return j

    j.image(image_path('vep'))
    j.storage('50G')
    j.cpu(16)

    if not isinstance(vcf, hb.Resource):
        vcf = b.read_input(str(vcf))

    if out_format == 'vcf':
        j.declare_resource_group(
            output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
        )
        output = j.output['vcf.gz']
    else:
        output = j.output

    # gcsfuse works only with the root bucket, without prefix:
    base_bucket_name = reference_path('vep_mount').drive
    data_mount = f'/{base_bucket_name}'
    j.cloudfuse(base_bucket_name, str(data_mount))
    vep_dir = f'{data_mount}/vep/GRCh38'
    loftee_conf = {
        'loftee_path': '$LOFTEE_PLUGIN_PATH',
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
    }

    authenticate_cloud_credentials_in_job(j)
    cmd = f"""\
    ls {vep_dir}
    ls {vep_dir}/vep

    LOFTEE_PLUGIN_PATH=/root/micromamba/share/ensembl-vep-105.0-1
    FASTA={vep_dir}/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

    vep \\
    --format vcf \\
    --{out_format} {'--compress_output bgzip' if out_format == 'vcf' else ''} \\
    -o {output} \\
    -i {vcf} \\
    --everything \\
    --allele_number \\
    --minimal \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache {vep_dir}/vep/ \\
    --dir_plugins $LOFTEE_PLUGIN_PATH \\
    --fasta $FASTA \\
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())}
    """
    if out_format == 'vcf':
        cmd += f'tabix -p vcf {output} '

    j.command(cmd)
    if out_path:
        b.write_output(j.output, str(out_path).replace('.vcf.gz', ''))
    return j
