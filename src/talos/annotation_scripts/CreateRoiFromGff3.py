#!/usr/bin/env python3

"""
Parses a GFF3 file, and generates a BED file of gene regions, plus padding, generating the following columns:
- chrom
- start
- end
- gene details
"""

import gzip
import re
from argparse import ArgumentParser


CHROM_INDEX = 0
RESOURCE_INDEX = 1
TYPE_INDEX = 2
START_INDEX = 3
END_INDEX = 4
DETAILS_INDEX = 8

# +/- this is added to each gene region, this default can be overridden
FLANKING_REGION = 2000

# regular expressions to parse out sections of the GFF3 annotations
GENE_ID_RE = re.compile(r'gene:(ENSG\d+);')
GENE_NAME_RE = re.compile(r'Name=([\w-]+);')

TYPES_TO_KEEP: set[str] = {'gene', 'ncRNA_gene', 'snRNA'}


def main(gff3_file: str, unmerged_output: str, merged_output: str, flanking: int = 2000):
    """
    Read the GFF3 file, and generate a BED file of gene regions, plus padding
    Args:
        gff3_file (str): path to the GFF3 file
        unmerged_output (str): path to the intended BED output file
        merged_output (str): path to generate a BED file with merged overlapping rows
        flanking (int): number of bases to add before and after each gene
    """

    unmerged_lines = generate_bed_lines(gff3_file, unmerged_output, flanking)
    merge_output(unmerged_lines, merged_output)


def generate_bed_lines(
    gff3_file: str,
    output: str,
    flanking: int = FLANKING_REGION,
) -> list[tuple[str, int, int]]:
    """
    generate the new BED file, and return the lines as a list of lists for merging
    """
    output_lines: list[tuple[str, int, int]] = []
    # open and iterate over the GFF3 file
    with gzip.open(gff3_file, 'rt') as handle, open(output, 'w') as write_handle:
        for line in handle:
            # skip over headers and dividing lines
            if line.startswith('#'):
                continue

            line_as_list = line.rstrip().split('\t')

            # skip over non-genes (e.g. pseudogenes, ncRNA)
            # only focus on Ensembl genes/transcripts
            if line_as_list[TYPE_INDEX] not in TYPES_TO_KEEP or 'ensembl' not in line_as_list[RESOURCE_INDEX]:
                continue

            # extract the gene name from the details field
            # allowing for some situations that don't work,
            # e.g. ENSG00000225931 (novel transcript, to be experimentally confirmed)
            # search for ID and transcript separately, ordering not guaranteed
            gene_name_match = GENE_NAME_RE.search(line_as_list[DETAILS_INDEX])
            gene_id_match = GENE_ID_RE.search(line_as_list[DETAILS_INDEX])
            if gene_id_match and gene_name_match:
                gene_name = gene_name_match.group(1)
                gene_id = gene_id_match.group(1)
            else:
                print(f'Failed to extract gene name from {line_as_list[DETAILS_INDEX]}')
                continue

            # write the line to the output
            output_list = [
                f'chr{line_as_list[CHROM_INDEX]}',
                str(int(line_as_list[START_INDEX]) - flanking),
                str(int(line_as_list[END_INDEX]) + flanking),
                f'{gene_id};{gene_name}',
            ]
            write_handle.write('\t'.join(output_list) + '\n')
            output_lines.append(
                (
                    f'chr{line_as_list[CHROM_INDEX]}',
                    int(line_as_list[START_INDEX]) - flanking,
                    int(line_as_list[END_INDEX]) + flanking,
                ),
            )
    return output_lines


def merge_output(
    unmerged_lines: list[tuple[str, int, int]],
    output: str,
    flanking: int = FLANKING_REGION,
):
    """
    take the lines, resolve overlapping regions, and write out to a new file

    Args:
        unmerged_lines (list): input of region-per-gene lines
        output (str): filepath to generate
        flanking (int): merge regions within this distance

    Returns:

    """
    contig = None
    start = None
    end = None
    with open(output, 'w') as handle:
        for this_chrom, this_start, this_end in unmerged_lines:
            if contig is None:
                contig = this_chrom
                start = this_start
                end = this_end
                continue

            # change in contig = write the previous block, and reset values
            if this_chrom != contig:
                handle.write(f'{contig}\t{start}\t{end}\n')
                contig = this_chrom
                start = this_start
                end = this_end
                continue

            # adjacent blocks close enough, merge them, but don't write
            if this_start - end < flanking:
                end = this_end

            # far enough apart, write the previous block and reset
            else:
                handle.write(f'{contig}\t{start}\t{end}\n')
                start = this_start
                end = this_end

        # and write the final line
        handle.write(f'{contig}\t{start}\t{end}\n')


def cli_main():
    parser = ArgumentParser()
    parser.add_argument(
        '--gff3',
        help='Path to the compressed GFF3 file',
        required=True,
    )
    parser.add_argument(
        '--unmerged_output',
        help='Path to output file, one line per gene with IDs',
        required=True,
    )
    parser.add_argument(
        '--merged_output',
        help='Path to output file, regions merged',
        required=False,
    )
    parser.add_argument(
        '--flanking',
        help='Regions to add to each gene',
        default=FLANKING_REGION,
    )
    args = parser.parse_args()
    main(
        gff3_file=args.gff3,
        unmerged_output=args.unmerged_output,
        merged_output=args.merged_output,
        flanking=args.flanking,
    )


if __name__ == '__main__':
    cli_main()
