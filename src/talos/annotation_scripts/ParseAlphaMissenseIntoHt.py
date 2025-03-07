#! /usr/bin/env python3

"""
Expects the from-source AlphaMissense table, available from:
- https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_hg38.tsv.gz
or
- https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz

Script for parsing the compressed tsv file of AlphaMissense results into a Hail Table
This Hail Table can be used for integrating AM scores into VCFs annotated with VEP or BCFtools

The AlphaMissense data is publicly available, and this is just a parser for that TSV

Assumption: the first lines of the AM data are still
# Copyright 2023 DeepMind Technologies Limited
#
# Licensed under CC BY-NC-SA 4.0 license

Process:
1. read through the compressed data, skip non-pathogenic entries
2. write the pathogenic entries back out to a reduced schema
3. parse that data as a Hail Table using specific variant types
4. write the Hail Table
"""

import gzip
import json
from argparse import ArgumentParser

import hail as hl


def process_header(final_header_line: str) -> dict[str, int]:
    """
    the TSV format is a little different in the all-isoforms vs. main transcript only
    this method determines the column indexes to use based on the header content
    we could determine preset columns by filename/flag, but... that's more burden on operator

    Args:
        final_header_line (str): the header line from the TSV file
    """
    # remove newline and hash, lowercase, split into a list
    broken_line = final_header_line.rstrip().replace('#', '').lower().split()

    return {
        'chrom': broken_line.index('chrom'),
        'pos': broken_line.index('pos'),
        'ref': broken_line.index('ref'),
        'alt': broken_line.index('alt'),
        'transcript': broken_line.index('transcript_id'),
        'am_pathogenicity': broken_line.index('am_pathogenicity'),
        'am_class': broken_line.index('am_class'),
    }


def filter_for_pathogenic_am(input_file: str, intermediate_file: str):
    """
    read the tsv file, skim for pathogenic entries, then write out to a new file

    Args:
        input_file ():
        intermediate_file ():
    """

    headers = ['chrom', 'pos', 'ref', 'alt', 'transcript', 'am_pathogenicity', 'am_class']

    # empty dictionary to contain the target indexes
    header_indexes: dict[str, int] = {}
    with gzip.open(input_file, 'rt') as read_handle, open(intermediate_file, 'w') as write_handle:
        for line in read_handle:
            # skip over the headers
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    # set the indexes (isoform and main-only have different columns)
                    header_indexes = process_header(line)
                continue

            if not header_indexes:
                raise ValueError('No header line was identified, columns are a mystery')

            # skip over everything except pathogenic
            if 'pathogenic' not in line:
                continue

            content = line.rstrip().split()

            # grab all the content
            content_dict: dict[str, str | float] = {key: content[header_indexes[key]] for key in headers}

            # trim transcripts
            content_dict['transcript'] = str(content_dict['transcript']).split('.')[0]

            # convert the AM score to a float, and pos to an int
            content_dict['pos'] = int(content_dict['pos'])
            content_dict['am_pathogenicity'] = float(content_dict['am_pathogenicity'])
            write_handle.write(f'{json.dumps(content_dict)}\n')


def json_to_hail_table(json_file: str, new_ht: str):
    """
    take a previously created JSON file and ingest it as a Hail Table
    requires an initiated Hail context

    Args:
        json_file ():
        new_ht ():
    """

    # define the schema for each written line
    schema = hl.dtype(
        'struct{chrom:str,pos:int32,ref:str,alt:str,transcript:str,am_pathogenicity:float64,am_class:str}',
    )

    # import as a hail table, force=True as this isn't Block-Zipped so all read on one core
    # We also provide a full attribute schema
    ht = hl.import_table(json_file, types={'f0': schema}, force=True, no_header=True)
    ht = ht.transmute(**ht.f0)

    # combine the two alleles into a single list
    ht = ht.transmute(locus=hl.locus(contig=ht.chrom, pos=ht.pos), alleles=[ht.ref, ht.alt])
    ht = ht.key_by('locus', 'alleles')
    ht.write(new_ht)
    ht.describe()


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--am_tsv', help='path to the AM tsv.gz file')
    parser.add_argument('--ht_out', help='path to write a new Hail Table')
    args, unknown = parser.parse_known_args()

    if unknown:
        raise ValueError(unknown)
    main(alpha_m_file=args.am_tsv, ht_path=args.ht_out)


def main(alpha_m_file: str, ht_path: str):
    """
    takes the path to an AlphaMissense TSV, reorganises it into a Hail Table

    Args:
        alpha_m_file ():
        ht_path ():
    """

    # generate a random file name so that we don't overwrite anything consistently
    random_intermediate_file: str = 'temp.json'

    # generate a new tsv of just pathogenic entries
    filter_for_pathogenic_am(alpha_m_file, random_intermediate_file)

    hl.default_reference('GRCh38')

    # now ingest as HT and re-jig some fields
    json_to_hail_table(random_intermediate_file, ht_path)


if __name__ == '__main__':
    cli_main()
