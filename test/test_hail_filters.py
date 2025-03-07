"""
unit testing collection for the hail MT methods
"""

import hail as hl
import pytest

from talos.RunHailFiltering import filter_matrix_by_ac, filter_on_quality_flags


@pytest.mark.parametrize(  # needs clinvar
    'ac,an,clinvar,threshold,rows',
    [
        ([1], 1, 0, 0.01, 1),
        ([6], 1, 0, 0.01, 0),
        ([6], 1, 1, 0.01, 1),
        ([6], 70, 0, 0.1, 1),
        ([50], 999999, 0, 0.01, 1),
        ([50], 50, 0, 0.01, 0),
        ([50], 50, 1, 0.01, 1),
    ],
)
def test_ac_filter_no_filt(
    ac: int,
    an: int,
    clinvar: int,
    threshold: float,
    rows: int,
    make_a_mt: hl.MatrixTable,
):
    """
    run tests on the ac filtering method
    check that a clinvar pathogenic overrides the AC test
    """
    matrix = make_a_mt.annotate_rows(
        info=make_a_mt.info.annotate(
            clinvar_talos=clinvar,
            AC=ac,
            AN=an,
        ),
    )

    assert filter_matrix_by_ac(matrix, threshold).count_rows() == rows


@pytest.mark.parametrize(
    'filters,clinvar,length',
    [
        (hl.empty_set(hl.tstr), 0, 1),
        (hl.literal({'fail'}), 0, 0),
        (hl.literal({'fail'}), 1, 1),
        (hl.literal({'VQSR'}), 0, 0),
        (hl.literal({'VQSR'}), 1, 1),
    ],
)
def test_filter_on_quality_flags(
    filters: hl.set,
    clinvar: int,
    length: int,
    make_a_mt: hl.MatrixTable,
):
    """
    annotate filters and run tests
    """
    # to add new alleles, we need to scrub alleles from the key fields
    anno_matrix = make_a_mt.key_rows_by('locus')
    anno_matrix = anno_matrix.annotate_rows(
        filters=filters,
        info=anno_matrix.info.annotate(clinvar_talos=clinvar),
    )
    assert filter_on_quality_flags(anno_matrix).count_rows() == length
