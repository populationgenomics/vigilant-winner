"""
unit testing collection for the hail MT methods
"""

import hail as hl
import pytest

from talos.RunHailFiltering import filter_matrix_by_ac, filter_on_quality_flags, filter_to_well_normalised


@pytest.mark.parametrize(  # needs clinvar
    'ac,an,threshold,rows',
    [
        ([1], 1, 0.01, 1),
        ([6], 1, 0.01, 0),
        ([6], 70, 0.1, 1),
        ([50], 999999, 0.01, 1),
        ([50], 50, 0.01, 0),
    ],
)
def test_ac_filter_no_filt(
    ac: int,
    an: int,
    threshold: float,
    rows: int,
    make_a_mt: hl.MatrixTable,
):
    """
    run tests on the ac filtering method
    """
    matrix = make_a_mt.annotate_rows(info=make_a_mt.info.annotate(AC=ac, AN=an))
    assert filter_matrix_by_ac(matrix, threshold).count_rows() == rows


@pytest.mark.parametrize(
    'filters,clinvar,svdb,exomiser,length',
    [
        (hl.empty_set(hl.tstr), 0, 0, 'missing', 1),
        (hl.literal({'fail'}), 0, 0, 'missing', 0),
        (hl.literal({'fail'}), 0, 1, 'missing', 1),
        (hl.literal({'fail'}), 0, 0, 'present', 1),
        (hl.literal({'fail'}), 1, 0, 'missing', 1),
        (hl.literal({'VQSR'}), 0, 0, 'missing', 0),
        (hl.literal({'VQSR'}), 1, 0, 'missing', 1),
    ],
)
def test_filter_on_quality_flags(
    filters: hl.set,
    clinvar: int,
    svdb: str,
    exomiser: str,
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
        info=anno_matrix.info.annotate(
            clinvar_talos=clinvar,
            categorybooleansvdb=svdb,
            categorydetailsexomiser=exomiser,
        ),
    )
    assert filter_on_quality_flags(anno_matrix).count_rows() == length


@pytest.mark.parametrize(
    'alleles,length',
    [
        (hl.literal(['A', 'C']), 1),
        (hl.literal(['A', '*']), 0),
        (hl.literal(['A', 'C', 'G']), 0),
    ],
)
def test_filter_to_well_normalised(alleles, length, make_a_mt):
    """
    checks the allele-level tests
    """
    # to add new alleles, we need to scrub alleles from the key fields
    anno_matrix = make_a_mt.key_rows_by('locus')
    anno_matrix = anno_matrix.annotate_rows(alleles=alleles)
    assert filter_to_well_normalised(anno_matrix).count_rows() == length
