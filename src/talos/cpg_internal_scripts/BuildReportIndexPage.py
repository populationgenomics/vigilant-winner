"""
track down the latest version of all reports
generate an index HTML page with links to all reports
"""

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import jinja2
from cloudpathlib.anypath import to_anypath
from metamist.graphql import gql, query

from talos.static_values import get_logger


DATE_REGEX = re.compile(r'(\d{4}-\d{2}-\d{2})')

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'
PROJECT_QUERY = gql(
    """
    query MyQuery {
        myProjects {
            dataset
        }
    }
    """,
)
REPORT_QUERY = gql(
    """
    query MyQuery($project: String!) {
        project(name: $project) {
            analyses(active: {eq: true}, type:  {eq: "aip-report"}) {
                outputs
                meta
                timestampCompleted
            }
        }
    }
    """,
)

WEB_BASE = 'gs://cpg-{}-main-web'
WEB_URL_BASE = 'https://main-web.populationgenomics.org.au/{}'
INDEX_HOME = 'gs://cpg-common-test-web/reanalysis/{}'


script_logger = get_logger(logger_name=__file__)


@dataclass
class Report:
    """
    generic object for storing report details
    """

    dataset: str
    address: str
    genome_or_exome: str
    date: str
    title: str


@lru_cache(1)
def get_my_projects() -> set[str]:
    """
    queries metamist for projects I have access to,
    returns the dataset names
    """
    response: dict[str, Any] = query(PROJECT_QUERY)
    all_projects = {dataset['dataset'] for dataset in response['myProjects']}
    script_logger.info(f'Running for projects: {", ".join(sorted(all_projects))}')
    return all_projects


def get_project_analyses(project: str) -> dict[str, str]:
    """
    find all the active analysis entries for this project
    we only want one regular report, the latest

    Args:
        project (str): project to query for
    """

    # we no longer generate 'latest' reports - HTML limitations have been removed
    project_reports: dict[str, str] = {'exome': '', 'genome': ''}
    all_analyses = query(REPORT_QUERY, variables={'project': project})['project']['analyses']
    for analysis in all_analyses:
        # get the type or skip (outdated)
        if not (st := analysis['meta'].get('sequencing_type')):
            continue

        # get the output path, allow for old analysis entries
        output_path = analysis['outputs'] if isinstance(analysis['outputs'], str) else analysis['outputs']['path']

        project_reports[st] = output_path

    return project_reports


def main() -> None:
    """
    finds all existing reports, generates an HTML file
    """

    parsed_reports = {cohort: get_project_analyses(cohort) for cohort in get_my_projects()}

    report_list: list[Report] = []

    for cohort, cohort_results in parsed_reports.items():
        for sequencing_type, report_path in cohort_results.items():
            # general - only one of these
            if report_path:
                this_file_name = Path(report_path).name
                trimmed_path = report_path.rstrip(this_file_name).rstrip('/')

                for entry in list(map(str, to_anypath(trimmed_path).glob('*.html'))):
                    report_address = entry.replace(WEB_BASE.format(cohort), WEB_URL_BASE.format(cohort))
                    report_name = entry.split('/')[-1]
                    if report_date := DATE_REGEX.search(report_address):
                        report_list.append(
                            Report(
                                dataset=cohort,
                                address=report_address,
                                genome_or_exome=sequencing_type,
                                date=report_date.group(1),
                                title=report_name,
                            ),
                        )

    html_from_reports(report_list, 'aip_index.html')


def html_from_reports(reports: list[Report], title: str):
    """
    build some HTML
    Args:
        reports (list[Report]): list of reports to build HTML for
        title (str): title of the page
    """

    # smoosh into a list for the report context - all reports sortable by date
    template_context = {'reports': reports}

    # build some HTML
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('report_index.html.jinja')
    content = template.render(**template_context)

    # write to common web bucket - either attached to a single dataset, or communal
    write_index_to = to_anypath(INDEX_HOME.format(title))
    get_logger().info(f'Writing {title} to {write_index_to}')
    write_index_to.write_text('\n'.join(line for line in content.split('\n') if line.strip()))


if __name__ == '__main__':
    get_logger(__file__).info('Fetching all reports')
    main()
