"""
Takes an input MT, and extracts a VCF-format representation.

This is currently required as the end-to-end CPG pipeline doesn't currently
store intermediate files. To simulate workflows running on VCF files, we
have to regenerate a VCF representation from a MT.

Hard coded additional header file for VQSR content
When Hail extracts a VCF from a MT, it doesn't contain any custom field
definitions, e.g. 'VQSR' as a Filter field. This argument allows us to
specify additional lines which are required to make the final output valid
within the VCF specification

If the file was not processed with VQSR, there are no negatives to including
this additional header line

New behaviour - by default region filter prior to VCf export
"""

import gzip
import logging
import sys
from argparse import ArgumentParser

import requests
import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch, output_path


CANON_CHROMS = [f'chr{chrom}' for chrom in range(1, 23)] + ['X', 'Y', 'M']
# path for downloading GenCode GTF file
GENCODE_GTF_URL = (
    'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/'
    'release_{gencode_release}/gencode.v{gencode_release}.annotation.gtf.gz'
)

# distance between adjacent genes before they're merged into a single span
# compromise here between reducing regions and reducing intergenic capture
# n == 1000: 26400 BED rows
# n == 5000: 17500 BED rows
BED_MERGE_MARGIN = 5000
LOCAL_BED = 'localfile.bed'
LOCAL_GTF = 'localfile.gtf.gz'


def download_gencode(gencode_release: str = '44'):
    """
    Download the GTF file from GENCODE
    Args:
        gencode_release (str): Which gencode release do you want?
    Returns:
        str - path to localised GTF file
    """
    gtf_path = GENCODE_GTF_URL.format(gencode_release=gencode_release)
    gz_stream = requests.get(gtf_path, stream=True)
    with open(LOCAL_GTF, 'wb') as f:
        f.writelines(gz_stream)
    gz_stream.close()


def parse_gtf_from_local(gtf_override: str | None = None):
    """
    Read over the localised GTF and parse into a BED file
    This is done by looping over all gene entries, and
    condensing overlapping adjacent gene region definitions

    Args:
        gtf_override (str|None): non-default GTF to use
    Returns:
        the path to an interval merged BED file
    """

    gtf_file = gtf_override or LOCAL_GTF

    logging.info(f'Loading {gtf_file}')

    def strip_from_list(val_list: list) -> tuple[str, int, int]:
        """
        quick method to reduce line count
        """
        chrom_string = val_list[0]
        start_int = max(int(val_list[3]) - BED_MERGE_MARGIN, 1)
        end_int = min(
            int(val_list[4]) + BED_MERGE_MARGIN,
            hl.get_reference('GRCh38').lengths[chrom_string],
        )
        return chrom_string, start_int, end_int

    out_rows = []
    cur_chr = None
    start = None
    end = None
    with gzip.open(gtf_file, 'rt') as gencode_file:

        # iterate over this file and do all the things
        for i, line in enumerate(gencode_file):
            line = line.rstrip('\r\n')
            if not line or line.startswith('#') or ('gene' not in line):
                continue

            # example line of interest
            # chr1    HAVANA  gene    11869   14409   .       +
            fields = line.split('\t')
            if fields[2] != 'gene':
                continue

            # first line
            if cur_chr is None:
                cur_chr, start, end = strip_from_list(fields)
                continue

            elif cur_chr != fields[0]:
                out_rows.append([cur_chr, start, end])
                cur_chr, start, end = strip_from_list(fields)
                continue

            this_chr, this_start, this_end = strip_from_list(fields)
            if this_start < end or this_end < end:
                start = min(this_start, start)
                end = max(this_end, end)
            else:
                out_rows.append([cur_chr, start, end])
                cur_chr, start, end = this_chr, this_start, this_end

        out_rows.append([cur_chr, start, end])

    # save as a BED file, return path to the BED file
    with open(LOCAL_BED, 'w', encoding='utf-8') as bed_file:
        for row in out_rows:
            bed_file.write('\t'.join([str(x) for x in row]) + '\n')


def main(mt_path: str, write_path: str):
    """
    takes an input MT, and reads it out as a VCF
    inserted new conditions to minimise the data produced

    Args:
        mt_path ():
        write_path ():
    """
    init_batch()

    mt = hl.read_matrix_table(mt_path)

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    additional_cloud_path = output_path('additional_header.txt', 'tmp')

    with to_path(additional_cloud_path).open('w') as handle:
        handle.write('##FILTER=<ID=VQSR,Description="VQSR triggered">')

    # remove potentially problematic field from gVCF
    if 'gvcf_info' in mt.row_value:
        mt = mt.drop('gvcf_info')

    # apply region filtering:
    # download the GTF file
    # parse gene regions
    # collapse overlapping regions
    # write as a BED file
    # read in and annotate the MT
    # filter by defined intervals
    download_gencode('44')
    parse_gtf_from_local()
    interval_table = hl.import_bed(LOCAL_BED, reference_genome='GRCh38')
    filtered_mt = mt.filter_rows(hl.is_defined(interval_table[mt.locus]))

    # filter out non-variant rows
    filtered_mt = hl.variant_qc(filtered_mt)
    filtered_mt = filtered_mt.filter_rows(filtered_mt.variant_qc.n_non_ref > 0)

    # required to repeat the MT -> VCF -> MT cycle
    if 'info' in filtered_mt.row_value and 'AC' in filtered_mt.info:
        if isinstance(filtered_mt['info']['AC'], hl.Int32Expression):
            filtered_mt = filtered_mt.annotate_rows(
                info=filtered_mt.info.annotate(AC=[filtered_mt.info.AC])
            )
        else:
            filtered_mt = filtered_mt.annotate_rows(
                info=filtered_mt.info.annotate(AC=[1])
            )

    hl.export_vcf(
        filtered_mt,
        write_path,
        append_to_header=additional_cloud_path,
        tabix=True,
    )


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    parser = ArgumentParser()
    parser.add_argument('--input', type=str, help='input MatrixTable path')
    parser.add_argument('--output', type=str, help='path to write VCF out to')
    args = parser.parse_args()
    main(mt_path=args.input, write_path=args.output)
