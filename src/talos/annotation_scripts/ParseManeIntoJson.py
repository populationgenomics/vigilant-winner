#!/usr/bin/env python3

"""
Reformat the MANE summary file into a Hail Table
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/README.txt
"""

import gzip
import json
from argparse import ArgumentParser
from csv import DictReader
from collections import defaultdict

import hail as hl


def mane_to_json(input_path: str, output_path: str):
    """
    Reformat the MANE summary file into a JSON file
    I can't quite grok the annotation from table, so I'm doing the dumb version
    """
    transcript_dict: dict[str, dict[str, str]] = defaultdict(dict)

    with gzip.open(input_path, 'rt') as handle:
        reader = DictReader(handle, delimiter='\t')
        for line in reader:
            enst = line['Ensembl_nuc'].split('.')[0]  # ENST
            ensg = line['Ensembl_Gene'].split('.')[0]  # ENSG
            ensp = line['Ensembl_prot'].split('.')[0]  # ENSP
            symbol = line['symbol'].split('.')[0]  # symbol
            nm_id = line['RefSeq_nuc'].split('.')[0]  # NMID
            mane = line['MANE_status']  # "MANE Select" or "MANE Plus Clinical"

            transcript_dict[enst] = {
                'mane': mane,
                'ensp': ensp,
                'nm_id': nm_id,
                'ensg': ensg,
                'symbol': symbol,
            }

    with open(output_path, 'w') as handle:
        json.dump(transcript_dict, handle, indent=2)


def mane_to_ht(input_path: str, output_path: str):
    """
    Reformat the MANE summary file into a Hail Table
    """

    temp_tsv = 'temp.tsv'
    # read the file
    with gzip.open(input_path, 'rt') as handle, open(temp_tsv, 'w') as write_handle:
        write_handle.write('\t'.join(['enst', 'ensg', 'ensp', 'symbol', 'nm_id', 'mane']) + '\n')
        reader = DictReader(handle, delimiter='\t')
        for line in reader:
            line_elements = [
                line['Ensembl_nuc'].split('.')[0],  # ENST
                line['Ensembl_Gene'].split('.')[0],  # ENSG
                line['Ensembl_prot'].split('.')[0],  # ENSP
                line['symbol'].split('.')[0],  # symbol
                line['RefSeq_nuc'].split('.')[0],  # NMID
                line['MANE_status'],  # "MANE Select" or "MANE Plus Clinical"
            ]
            write_handle.write('\t'.join(line_elements) + '\n')

    # import as a hail table, force=True as this isn't Block-Zipped so all read on one core
    # No data types provided, all string columns
    ht = hl.import_table(temp_tsv, force=True)
    ht = ht.key_by('enst')
    ht.write(output_path, overwrite=True)
    ht.describe()
    ht.show()


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the MANE summary file')
    parser.add_argument('--output', help='Path to write the resulting HT')
    parser.add_argument('--format', choices=['json', 'ht'], default='json')
    args = parser.parse_args()

    if args.format == 'json':
        mane_to_json(input_path=args.input, output_path=args.output)
    elif args.format == 'ht':
        mane_to_ht(input_path=args.input, output_path=args.output)
    else:
        raise ValueError(f'Unrecognised format {args.format}')


if __name__ == '__main__':
    cli_main()
