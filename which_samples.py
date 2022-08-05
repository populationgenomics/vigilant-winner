"""
quick script to find the samples from a mt
"""
import hail as hl
from cpg_utils.hail_batch import init_batch, output_path
from sample_metadata.apis import SampleApi

if __name__ == '__main__':
    init_batch()
    mt = hl.read_matrix_table('gs://cpg-nagim-main/mt/v1-2.mt')
    tg_samples = (
        SampleApi().get_all_sample_id_map_by_internal('thousand-genomes').keys()
    )
    samples = set(mt.s.collect()).intersection(tg_samples)
    print(f'{len(samples)} will be kept')
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))
    hl.experimental.densify(mt)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref > 0)
    mt = mt.drop('variant_qc')
    mt.write(output_path('thousand_genomes_subset.mt'))
