"""
quick script to find the samples from a mt
"""
import hail as hl
from cpg_utils.hail_batch import init_batch


if __name__ == '__main__':
    init_batch()
    mt = hl.read_matrix_table('gs://cpg-nagim-main/mt/v1-2.mt')
    print(list(mt.s.collect()))
