"""
a script to post-process the Talos output data

produces a minimised representation of the output data,
containing only the data required for the SEQR app.

- Metadata, containing the description for each flag
- Results, containing per-individual results
    - Individual ID
        - Variant ID
            - Categories (list)
            - Support Variants (list)

Also produce a second version of the same, limited to phenotype-matches
"""

import json
from argparse import ArgumentParser

from talos.models import MiniForSeqr, MiniVariant, ResultData
from talos.static_values import get_logger


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='the input file to process')
    parser.add_argument('--output', help='the output file to write to')
    parser.add_argument('--pheno', help='the output file for phenotype-matched data', default=None)
    parser.add_argument('--external_map', help='mapping of internal to external IDs for seqr', default=None)
    args = parser.parse_args()

    main(input_file=args.input, output=args.output, ext_map=args.external_map)
    if args.pheno:
        main(input_file=args.input, output=args.pheno, ext_map=args.external_map, pheno_match=True)


def main(input_file: str, output: str, ext_map: str | None = None, pheno_match: bool = False):
    """
    reads in the input file, shrinks it, and writes the output file

    Args:
        input_file (str):
        output (str):
        ext_map (str): optional mapping of internal to external IDs for seqr
        pheno_match (bool): whether to limit to phenotype-matching variants
    """

    if pheno_match:
        get_logger().info('Limiting to phenotype-matching variants')

    with open(input_file, encoding='utf-8') as f:
        data = ResultData.model_validate(json.load(f))

    lil_data = MiniForSeqr(metadata={'categories': data.metadata.categories})
    ext_map_dict = None
    if ext_map:
        with open(ext_map, encoding='utf-8') as f:
            ext_map_dict = json.load(f)

    for individual, details in data.results.items():
        # optionally update to point to Seqr identities
        indi_id = ext_map_dict.get(individual, individual) if ext_map_dict else individual

        lil_data.results[indi_id] = {}
        for variant in details.variants:
            var_data = variant.var_data
            if pheno_match and not variant.panels.matched:
                continue
            lil_data.results[indi_id][var_data.info['seqr_link']] = MiniVariant(
                categories=variant.categories,
                support_vars=variant.support_vars,
            )

    if not any(lil_data.results.values()):
        get_logger().info('No results found')

    # write anyway, so as not to break the pipeline
    with open(output, 'w', encoding='utf-8') as f:
        f.write(MiniForSeqr.model_validate(lil_data).model_dump_json(indent=4))

    get_logger().info(f'Wrote output to {output}')


if __name__ == '__main__':
    cli_main()
