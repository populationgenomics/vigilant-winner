#!/usr/bin/env python3

"""
A standalone script to read some talos report files, and summarise the number of affected samples involved

Allows for the printing of affected participants grouped by a portion of the family ID
Examples are when the family  ID is prefixed with the year, e.g. 19DNAXXXX for 2019
The simple count function doesn't have this level of granularity

If there is no common prefix, we don't attempt this grouping
"""

import json

from argparse import ArgumentParser
from collections import defaultdict, Counter

from cpg_utils import to_path

from talos.models import ResultData

from talos.utils import read_json_from_path

MEAN_SLASH_SAMPLE = 'Mean/sample'


class NoVariantsFoundError(Exception):
    """raise if a report subset contains no data"""


def get_variant_summary(results: ResultData) -> dict:
    """
    Run the numbers across all variant categories
    Treat each primary-secondary comp-het pairing as one event
    i.e. the thing being counted here is the number of events
    which passed through the MOI process, not the absolute number
    of variants in the report

    Args:
        results (ResultData): the results object in full

    Returns:
        a dictionary summarising the categorised variants
    """

    # get the categories this report was aware of
    all_categories = results.metadata.categories.keys()

    ordered_categories = ['any', *all_categories]

    category_count: dict = {key: [] for key in ordered_categories}
    unique_variants: dict[str, set[str]] = {key: set() for key in ordered_categories}

    # this keeps a count of the total instances of each variant
    global_count: dict[str, int] = defaultdict(int)

    for sample_data in results.results.values():
        sample_variants: dict[str, set[str]] = {key: set() for key in ordered_categories}

        # iterate over the list of variants
        for variant in sample_data.variants:
            var_string = variant.var_data.coordinates.string_format
            global_count[var_string] += 1
            unique_variants['any'].add(var_string)
            sample_variants['any'].add(var_string)

            # find all categories associated with this variant
            # for each category, add to corresponding list and set
            for category_value in variant.categories:
                unique_variants[category_value].add(var_string)
                sample_variants[category_value].add(var_string)

        # update the global lists with per-sample counts
        for key, key_list in category_count.items():
            key_list.append(len(sample_variants[key]))

    summary_dicts = {
        key: {
            'Description': results.metadata.categories.get(key, 'All Variants'),
            'Total': sum(category_count[key]),
            'Unique': len(unique_variants[key]),
            'Peak #/sample': max(category_count[key]),
            MEAN_SLASH_SAMPLE: sum(category_count[key]) / len(category_count[key]),
        }
        for key in ordered_categories
    }

    # make a Counter object from the collected counts, identify the 10 most frequent
    summary_dicts['most_common'] = dict(Counter(global_count).most_common(10))

    # this can fail if there are no categorised variants... at all
    if not summary_dicts:
        raise NoVariantsFoundError('No categorised variants found')

    summary_dicts['samples_no_variants'] = category_count['any'].count(0)

    return summary_dicts


def main(input_path: str, output_path: str | None = None, prefix: int | None = None):
    """
    read the target report, and summarise the content:
      - the number of affected samples involved
      - the number of variants in each category
      - the number of samples with no variants

    Args:
        input_path (str): where to read the report from
        output_path (str): optional, if we want to write the output to a file
        prefix (int): optional, if we want to break down by family prefix
    """

    # read the report file, local or cloud
    report = read_json_from_path(input_path, return_model=ResultData)

    summarised_content: dict = {'family_breakdown': report.metadata.family_breakdown}

    if prefix:
        # set up a section in the dictionary for this
        summarised_content['family_breakdown']['grouped_by_prefix'] = defaultdict(int)
        for proband in report.results.values():
            summarised_content['family_breakdown']['grouped_by_prefix'][proband.metadata.ext_id[:prefix]] += 1

    summarised_content['variant_summary'] = get_variant_summary(report)

    print(json.dumps(summarised_content, indent=4))

    if output_path:
        # write the output to file
        with to_path(output_path).open('w') as handle:
            json.dump(summarised_content, handle, indent=4)


def cli_main():
    """
    CLI entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the input report file')
    parser.add_argument('--output', help='Where to write the output')
    parser.add_argument('--prefix', type=int, help='breakdown by family prefix', default=None)
    args = parser.parse_args()
    main(
        input_path=args.input,
        output_path=args.output,
        prefix=args.prefix,
    )


if __name__ == '__main__':
    cli_main()
