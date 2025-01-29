"""
Methods for taking the final output and generating static report content

This is another total rewrite, which tries to fit some resource-friendly
frontage onto the report, so that it loads in good time.

If there's a common prefix (e.g. by year), we split the data into sub-reports,
but we don't need to keep paring it down and down
"""

import re
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from os import makedirs
from os.path import join
from pathlib import Path
from typing import Any

import jinja2
import pandas as pd

from talos.config import config_retrieve
from talos.models import PanelApp, PanelDetail, ReportVariant, ResultData, SmallVariant, StructuralVariant
from talos.utils import get_logger, read_json_from_path

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'
MIN_REPORT_SIZE: int = 10
MAX_REPORT_SIZE: int = 200

# above this length we trim the actual bases to just an int
MAX_INDEL_LEN: int = 10

# regex pattern - number, number, not-number
KNOWN_YEAR_PREFIX = re.compile(r'\d{2}\D')
CDNA_SQUASH = re.compile(r'(?P<type>ins|del)(?P<bases>[ACGT]+)$')
MEAN_SLASH_SAMPLE = 'Mean/sample'
GNOMAD_SV_KEY = 'gnomad_v2.1_sv_svid'


def known_date_prefix_check(all_results: ResultData) -> list[str]:
    """
    Check for known date prefixes in the results

    Args:
        all_results (): the whole summary dataset

    Returns:
        a list of all found prefixes, or empty list
    """

    known_prefixes: dict[str, int] = defaultdict(int)
    for content in all_results.results.values():
        if match := KNOWN_YEAR_PREFIX.match(content.metadata.ext_id):
            known_prefixes[match.group()[0:2]] += 1
        else:
            get_logger().info('There is no consistent sample ID prefix')
            return []

    get_logger().info(f'Sample distribution by prefix: {dict(known_prefixes)}')
    return sorted(known_prefixes.keys())


def split_data_into_sub_reports(all_results: ResultData) -> list[tuple[ResultData, str, str]]:
    """
    Split the data into sub-reports, only if there's a common prefix (e.g. by year)
    Return a list of the ResultData subsets, output base path, and a subset identifier

    Args:
        all_results ():

    Returns:
        tuple: a list of tuples, each containing a ResultData object, the output base path, and a subset identifier
    """
    return_results: list[tuple[ResultData, str, str]] = []

    if prefixes := known_date_prefix_check(all_results):
        for prefix in prefixes:
            this_rd = ResultData(
                metadata=all_results.metadata,
                results={
                    sample: content
                    for sample, content in all_results.results.items()
                    if content.metadata.ext_id.startswith(prefix)
                },
                version=all_results.version,
            )
            get_logger().info(f'Found {len(this_rd.results)} with prefix {prefix}')
            return_results.append((this_rd, f'subset_{prefix}.html', prefix))
    return return_results


def cli_main():
    get_logger(__file__).info('Running HTML builder')
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to analysis results', required=True)
    parser.add_argument('--panelapp', help='PanelApp data', required=True)
    parser.add_argument('--output', help='Final HTML filename', required=True)
    parser.add_argument('--latest', help='Not in use')
    parser.add_argument('--split_samples', help='Not in use')
    args = parser.parse_args()

    if args.latest:
        get_logger(__file__).warning('"--latest" argument is not in use')

    if args.split_samples:
        get_logger(__file__).warning('"--split_samples" argument is not in use')

    main(results=args.input, panelapp=args.panelapp, output=args.output)


def main(results: str, panelapp: str, output: str):
    """

    Args:
        results (str): path to the MOI-tested results file
        panelapp (str): path to the panelapp data
        output (str): where to write the HTML file
    """

    report_output_dir = Path(output).parent

    results_object = read_json_from_path(results, return_model=ResultData)

    # we always make this main page - we need a reliable output path to generate analysis entries [CPG]
    html = HTMLBuilder(results_dict=results_object, panelapp_path=panelapp)
    # if this fails with a NoVariantsFoundException, there were no variants to present in the whole cohort
    # catch this, but fail gracefully so that the process overall is a success
    try:
        get_logger().debug(f'Writing whole-cohort categorised variants to {output}')
        # find the path to the output directory, and make an individual directory
        makedirs(join(report_output_dir, 'individuals'), exist_ok=True)
        html.write_html(output_filepath=output)
    except NoVariantsFoundError:
        get_logger().warning('No Categorised variants found in this whole cohort')

    # we only need to do sub-reports if we can delineate by year
    for data, report, prefix in split_data_into_sub_reports(results_object):
        html = HTMLBuilder(results_dict=data, panelapp_path=panelapp, subset_id=prefix)
        try:
            output_filepath = join(report_output_dir, report)
            get_logger().debug(f'Attempting to create {report} at {output_filepath}')
            makedirs(join(report_output_dir, f'individuals_{prefix}'), exist_ok=True)
            html.write_html(output_filepath=output_filepath)
        except NoVariantsFoundError:
            get_logger().info('No variants in that report, skipping')


class NoVariantsFoundError(Exception):
    """raise if a report subset contains no data"""


@dataclass
class DataTable:
    """
    Representation of a DataTables table that the Jinja2 templating system renders.
    """

    id: str
    columns: list[str]
    rows: list[Any]
    heading: str = ''
    description: str = ''


def variant_in_forbidden_gene(variant_obj: ReportVariant, forbidden_genes):
    """
    Check if gene id or gene symbol is on forbidden gene list
    """
    for gene_id in variant_obj.gene.split(','):
        if gene_id in forbidden_genes:
            return True

    if isinstance(variant_obj.var_data, StructuralVariant):
        return False

    # Allow for exclusion by Symbol too
    return any(tx_con['symbol'] in forbidden_genes for tx_con in variant_obj.var_data.transcript_consequences)


class HTMLBuilder:
    """
    Takes the input, makes the output
    """

    def __init__(self, results_dict: ResultData, panelapp_path: str, subset_id: str | None = None):
        """
        Args:
            results_dict (ResultData): the results object
            panelapp_path (str): where to read panelapp data from
            subset_id (str, optional): the subset ID to use for this report
        """
        # ID to use if this is a subset report
        self.subset_id = subset_id

        # get a hold of the base panel ID we're using
        # this is used to differentiate between new in base and new in other
        self.base_panel: int = config_retrieve(['GeneratePanelData', 'default_panel'], 137)

        self.panelapp: PanelApp = read_json_from_path(panelapp_path, return_model=PanelApp)

        # If it exists, read the forbidden genes as a list
        self.forbidden_genes = config_retrieve(['GeneratePanelData', 'forbidden_genes'], [])
        assert isinstance(self.forbidden_genes, list)
        get_logger().warning(f'There are {len(self.forbidden_genes)} forbidden genes')

        # Use config to find CPG-to-Seqr ID JSON; allow to fail
        self.seqr: dict[str, str] = {}

        if seqr_path := config_retrieve(['CreateTalosHTML', 'seqr_lookup'], ''):
            self.seqr = read_json_from_path(seqr_path, default=self.seqr)
            assert isinstance(self.seqr, dict)

            # Force user to correct config file if seqr URL/project are missing
            for seqr_key in ['seqr_instance', 'seqr_project']:
                assert config_retrieve(['CreateTalosHTML', seqr_key], False), f'Seqr key absent: {seqr_key}'

        # Optionally read in the labels file
        # This file should be a nested dictionary of sample IDs and variant identifiers
        # with a list of corresponding label values, e.g.:
        # ruff: noqa: ERA001
        # {
        #     "sample1": {
        #         "1-123456-A-T": ["label1", "label2"],
        #         "1-123457-A-T": ["label1"]
        #     },
        # }
        self.ext_labels: dict[str, dict] = config_retrieve(['CreateTalosHTML', 'external_labels'], {})
        assert isinstance(self.ext_labels, dict)

        self.metadata = results_dict.metadata
        self.panel_names = {panel.name for panel in self.metadata.panels}

        # if a variant has any filters (AB Ratio, Low Depth), the variant will only be reported if it's in this list
        # needs to go here because it's used in the loop below
        self.allow_filters: set[str] = set(config_retrieve(['CreateTalosHTML', 'allow_filters'], []))

        # Process samples and variants
        self.samples: list[Sample] = []
        for sample, content in results_dict.results.items():
            # skip for now if there's nothing to show
            if not content.variants:
                continue

            self.samples.append(
                Sample(
                    name=sample,
                    metadata=content.metadata,
                    variants=content.variants,
                    ext_labels=self.ext_labels.get(sample, {}),
                    html_builder=self,
                ),
            )

        self.samples.sort(key=lambda x: x.ext_id)

    def get_summary_stats(self) -> tuple[pd.DataFrame, list[str], list[dict]]:
        """
        Run the numbers across all variant categories
        Treat each primary-secondary comp-het pairing as one event
        i.e. the thing being counted here is the number of events
        which passed through the MOI process, not the absolute number
        of variants in the report
        """
        ordered_categories = ['any', *list(config_retrieve('categories', {}).keys())]
        category_count: dict[str, list[int]] = {key: [] for key in ordered_categories}
        unique_variants: dict[str, set[str]] = {key: set() for key in ordered_categories}

        samples_with_no_variants: list[str] = []
        ext_label_map: dict = self.ext_labels.copy() if self.ext_labels else {}

        for sample in self.samples:
            if len(sample.variants) == 0:
                samples_with_no_variants.append(sample.ext_id)

            sample_variants: dict[str, set[str]] = {key: set() for key in ordered_categories}

            # iterate over the list of variants
            for variant in sample.variants:
                var_string = variant.var_data.coordinates.string_format
                unique_variants['any'].add(var_string)
                sample_variants['any'].add(var_string)

                # find all categories associated with this variant
                # for each category, add to corresponding list and set
                for category_value in variant.categories:
                    sample_variants[category_value].add(var_string)
                    unique_variants[category_value].add(var_string)

                # remove any external labels associated with this sample/variant.
                if sample.name in ext_label_map:
                    ext_label_map[sample.name].pop(var_string, None)

            # update the global lists with per-sample counts
            for key, key_list in category_count.items():
                key_list.append(len(sample_variants[key]))

        # Extract the list of unused ext labels
        unused_ext_labels = [
            {
                'sample': sample_id,
                'sample_ext': self.seqr.get(sample_id, sample_id),
                'variant': var_id,
                'labels': labels,
            }
            for sample_id, var_dict in ext_label_map.items()
            for var_id, labels in var_dict.items()
        ]

        summary_dicts = [
            {
                'Category': key,
                'Total': sum(category_count[key]),
                'Unique': len(unique_variants[key]),
                'Peak #/sample': max(category_count[key]),
                MEAN_SLASH_SAMPLE: sum(category_count[key]) / len(category_count[key]),
            }
            for key in ordered_categories
            if category_count[key]
        ]

        # this can fail if there are no categorised variants... at all
        if not summary_dicts:
            raise NoVariantsFoundError('No categorised variants found')

        my_df: pd.DataFrame = pd.DataFrame(summary_dicts)
        my_df[MEAN_SLASH_SAMPLE] = my_df[MEAN_SLASH_SAMPLE].round(3)

        # the table re-sorts when parsed into the DataTable
        # so this forced ordering doesn't work
        my_df.Category = my_df.Category.astype('category')
        my_df.Category = my_df.Category.cat.set_categories(ordered_categories)
        my_df = my_df.sort_values(by='Category')

        return my_df, samples_with_no_variants, unused_ext_labels

    def read_metadata(self) -> dict[str, pd.DataFrame]:
        """
        parses into a general table and a panel table
        """

        return {
            'Panels': pd.DataFrame(
                {'ID': panel.id, 'Version': panel.version, 'Name': panel.name} for panel in self.metadata.panels
            ),
            'Meta': pd.DataFrame(
                {
                    'Data': key.capitalize(),
                    'Value': self.metadata.__getattribute__(key),
                }
                for key in ['run_datetime', 'version']
            ),
            'Families': pd.DataFrame(
                [
                    {'family_size': fam_type, 'tally': fam_count}
                    for fam_type, fam_count in sorted(self.metadata.family_breakdown.items())
                ],
            ),
        }

    def write_html(self, output_filepath: str):
        """
        Uses the results to create the HTML tables
        writes all content to the output path

        Args:
            output_filepath (str): where to write the results to
        """

        # if no variants were found, this can fail with a NoVariantsFoundException error
        # we ignore that here, and catch it in the outer scope
        (summary_table, zero_cat_samples, unused_ext_labels) = self.get_summary_stats()

        template_context = {
            # 'metadata': self.metadata,
            'index_path': f'../{Path(output_filepath).name}',
            'run_datetime': self.metadata.run_datetime,
            'samples': self.samples,
            'seqr_url': config_retrieve(['CreateTalosHTML', 'seqr_instance'], ''),
            'seqr_project': config_retrieve(['CreateTalosHTML', 'seqr_project'], ''),
            # 'meta_tables': {},
            # 'forbidden_genes': sorted(self.forbidden_genes),
            # 'zero_categorised_samples': zero_cat_samples,
            # 'unused_ext_labels': unused_ext_labels,
            # 'summary_table': None,
            'report_title': 'Full Talos Report',
            # 'solved': self.solved,
            'type': 'whole_cohort',
        }

        # write all HTML content to the output file in one go
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
            autoescape=True,
            trim_blocks=True,
            keep_trailing_newline=False,
            lstrip_blocks=True,
        )
        template = env.get_template('index.html.jinja')
        content = template.render(**template_context)
        with open(output_filepath, 'w') as handle:
            handle.writelines(content)

        get_logger().info(f'Wrote {output_filepath}')

        outpath_name = Path(output_filepath).name

        for sample in template_context['samples']:
            if not sample.variants:
                continue

            report_address = output_filepath.replace(outpath_name, sample.report_url)

            get_logger().debug(f'Writing {report_address}')

            new_context = template_context | {
                'samples': [sample],
                'report_title': f'Talos Report for {sample.name}',
                'type': 'sample',
            }

            template = env.get_template('sample_index.html.jinja')
            content = template.render(**new_context)
            with open(report_address, 'w') as handle:
                handle.writelines(content)


class Sample:
    """
    Sample related logic
    """

    def __init__(
        self,
        name: str,
        metadata,
        variants: list[ReportVariant],
        ext_labels: dict[str, list[str]],
        html_builder: HTMLBuilder,
    ):
        indi_folder = f'individuals_{html_builder.subset_id}' if html_builder.subset_id is not None else 'individuals'
        self.metadata = metadata
        self.name = name
        self.family_id = metadata.family_id
        self.family_members = metadata.members
        self.phenotypes = metadata.phenotypes
        self.ext_id = metadata.ext_id
        self.panel_details = metadata.panel_details
        self.seqr_id = html_builder.seqr.get(name, None)
        self.report_url = f'{indi_folder}/{self.name}.html'

        # exclude any on the forbidden gene list
        variants = [var for var in variants if not variant_in_forbidden_gene(var, html_builder.forbidden_genes)]

        # drop variants which fail quality/depth filters, unless all filters are explicitly permitted
        variants = [var for var in variants if not (var.flags - html_builder.allow_filters)]

        self.variants = [
            Variant(
                report_variant,
                self,
                ext_labels.get(report_variant.var_data.coordinates.string_format, []),
                html_builder,
            )
            for report_variant in variants
            if not variant_in_forbidden_gene(report_variant, html_builder.forbidden_genes)
        ]

    def __str__(self):
        return self.name


class Variant:
    """
    Handle as much of per variant logic as we can here. Hopefully, this is just simple
    data munging and mapping operations.

    Try not to put presentation logic here - keep it in the jinja templates
    """

    def get_var_change(self) -> str:
        """
        Find the variant change for the variant
        - we want to truncate huge small variant InDels (ballooning column width)
           - e.g. LOLOLOLOLOLOLOLOLOLOLOLOLOLOLOLOLO->A -> del 34bp
        - SVs always need to be presented differently
           - e.g INS 4079bp
        """
        if isinstance(self.var_data, SmallVariant):
            if len(self.ref) > MAX_INDEL_LEN or len(self.alt) > MAX_INDEL_LEN:
                ref_len = len(self.ref)
                alt_len = len(self.alt)
                if ref_len > alt_len:
                    return f'del {ref_len - alt_len}bp'
                if ref_len == alt_len:
                    return f'complex delins {ref_len}bp'
                return f'ins {alt_len - ref_len}bp'

            return f'{self.ref}->{self.alt}'
        if isinstance(self.var_data, StructuralVariant):
            return f"{self.var_data.info['svtype']} {self.var_data.info['svlen']}bp"

        raise ValueError(f'Unknown variant type: {self.var_data.__class__.__name__}')

    def __init__(self, report_variant: ReportVariant, sample: Sample, ext_labels: list, html_builder: HTMLBuilder):
        self.var_data = report_variant.var_data
        self.var_type = report_variant.var_data.__class__.__name__
        self.chrom = report_variant.var_data.coordinates.chrom
        self.pos = report_variant.var_data.coordinates.pos
        self.ref = report_variant.var_data.coordinates.ref
        self.alt = report_variant.var_data.coordinates.alt
        self.change = self.get_var_change()
        self.categories = report_variant.categories
        self.first_tagged: str = report_variant.first_tagged
        self.support_vars = report_variant.support_vars
        self.warning_flags = report_variant.flags
        # these are the panel IDs which are matched based on HPO matching in PanelApp
        self.pheno_matches = {f'{name}({pid})' for pid, name in report_variant.panels.matched.items()}
        # these are the panel IDs we manually applied to this whole cohort
        self.forced_matches = {f'{name}({pid})' for pid, name in report_variant.panels.forced.items()}

        # collect all forced and matched panel IDs
        match_ids = set(report_variant.panels.forced.keys()).union(set(report_variant.panels.matched.keys())) - {
            html_builder.base_panel,
        }

        self.reasons = report_variant.reasons
        self.genotypes = report_variant.genotypes
        self.ext_labels = ext_labels
        # add the phenotype match date and HPO term id/labels
        self.phenotype_match_date = report_variant.date_of_phenotype_match
        self.phenotype_matches = report_variant.phenotype_labels

        # check if this variant is new in the base panel
        self.new_in_base_panel: bool = False

        # store if this variant is new in any of the other panels
        self.new_panels: set[str] = set()

        # List of (gene_id, symbol)
        self.genes: list[tuple[str, str]] = []
        for gene_id in report_variant.gene.split(','):
            gene_panelapp_entry = html_builder.panelapp.genes.get(gene_id, PanelDetail(symbol=gene_id))
            self.genes.append((gene_id, gene_panelapp_entry.symbol))

            # is this a new gene?
            new_panels = gene_panelapp_entry.new

            if html_builder.base_panel in new_panels:
                self.new_in_base_panel = True

            # now draw the rest of the owl
            self.new_panels.update(
                {f'{sample.metadata.panel_details[pid]}({pid})' for pid in new_panels.intersection(match_ids)},
            )

        # Summaries CSQ strings
        if isinstance(self.var_data, SmallVariant):
            (self.mane_csq, self.mane_hgvsps) = self.parse_csq()

        # pull up the highest AlphaMissense score, if present
        am_scores = (
            [
                float(csq['am_pathogenicity'])
                for csq in self.var_data.transcript_consequences
                if csq.get('am_pathogenicity')
            ]
            if isinstance(self.var_data, SmallVariant)
            else []
        )

        # make these accessible in the report/presentation
        self.var_data.info['alpha_missense_max'] = max(am_scores) if am_scores else 'missing'
        self.exomiser_results = report_variant.exomiser_results

        # this is the weird gnomad callset ID
        if (
            isinstance(self.var_data, StructuralVariant)
            and GNOMAD_SV_KEY in self.var_data.info
            and isinstance(self.var_data.info[GNOMAD_SV_KEY], str)
        ):
            self.var_data.info['gnomad_key'] = self.var_data.info[GNOMAD_SV_KEY].split('v2.1_')[-1]  # type: ignore[union-attr]

    def __str__(self) -> str:
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def parse_csq(self):
        """
        Parse CSQ variant string returning:
            - set of "consequences" from MANE transcripts
            - Set of variant effects in p. nomenclature (or c. if no p. is available)

        condense massive cdna annotations, e.g.
        c.4978-2_4978-1insAGGTAAGCTTAGAAATGAGAAAAGACATGCACTTTTCATGTTAATGAAGTGATCTGGCTTCTCTTTCTA
        """
        mane_consequences = set()
        mane_hgvsps = set()

        for csq in self.var_data.transcript_consequences:
            if 'consequence' not in csq:
                continue

            # if csq['mane_select'] or csq['mane_plus_clinical']:
            if csq['mane_select']:
                mane_consequences.update(csq['consequence'].split('&'))
                if csq['hgvsp']:
                    mane_hgvsps.add(csq['hgvsp'].split(':')[1])
                elif csq['hgvsc']:
                    hgvsc = csq['hgvsc'].split(':')[1]

                    # if massive indel base stretches are included, replace with a numerical length
                    if match := CDNA_SQUASH.search(hgvsc):
                        hgvsc.replace(match.group('bases'), str(len(match.group('bases'))))

                    mane_hgvsps.add(hgvsc)

        # simplify the consequence strings
        mane_consequences = ', '.join(_csq.replace('_variant', '').replace('_', ' ') for _csq in mane_consequences)
        mane_hgvsps = ', '.join(mane_hgvsps)

        return mane_consequences, mane_hgvsps


if __name__ == '__main__':
    cli_main()
