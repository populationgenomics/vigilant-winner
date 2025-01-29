"""
tests relating to the MOI filters
"""

from test.test_utils import TWO_EXPECTED
from unittest import mock

import pytest

from talos.models import Coordinates, SmallVariant
from talos.moi_tests import (
    BaseMoi,
    DominantAutosomal,
    MOIRunner,
    RecessiveAutosomalCH,
    RecessiveAutosomalHomo,
    XDominant,
    XPseudoDominantFemale,
    XRecessiveFemaleCH,
    XRecessiveFemaleHom,
    XRecessiveMale,
    check_for_second_hit,
)
from talos.utils import make_flexible_pedigree

TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
TEST_COORDS2 = Coordinates(chrom='2', pos=2, ref='G', alt='T')
TEST_COORDS_X_1 = Coordinates(chrom='X', pos=1, ref='G', alt='T')
TEST_COORDS_X_2 = Coordinates(chrom='X', pos=2, ref='G', alt='T')


@pytest.mark.parametrize(
    'first,comp_hets,sample,values',
    (
        ('', {}, '', []),  # no values
        ('', {}, 'a', []),  # sample not present
        ('', {'a': {'foo': []}}, 'a', []),  # var not present
        ('foo', {'a': {'foo': ['bar']}}, 'a', ['bar']),  # all values present
        ('foo', {'a': {'foo': ['bar', 'baz']}}, 'a', ['bar', 'baz']),  # all values present
    ),
)
def test_check_second_hit(first, comp_hets, sample, values):
    """
    quick test for the 2nd hit mechanic
    return all strings when the comp-het lookup contains:
        - the sample
        - the gene
        - the variant signature
    """

    assert check_for_second_hit(first_variant=first, comp_hets=comp_hets, sample=sample) == values


@pytest.mark.parametrize(
    'moi_string,filters',
    (
        ('Monoallelic', ['DominantAutosomal']),
        ('Mono_And_Biallelic', ['DominantAutosomal', 'RecessiveAutosomal']),
        ('Unknown', ['DominantAutosomal', 'RecessiveAutosomal']),
        ('Biallelic', ['RecessiveAutosomal']),
        ('Hemi_Mono_In_Female', ['XRecessive', 'XDominant']),
        ('Hemi_Bi_In_Female', ['XRecessive']),
    ),
)
def test_moi_runner(moi_string: str, filters: list[str], pedigree_path):
    """
    check that the right methods are associated with each MOI
    """
    test_runner = MOIRunner(pedigree=make_flexible_pedigree(pedigree_path), target_moi=moi_string)

    # the imported (uninstantiated) objects don't have __class__
    # and the instantiated objects don't have a __name__
    for filter1, filter2 in zip(test_runner.filter_list, filters, strict=True):
        assert filter2 in str(filter1.__class__)


def test_dominant_autosomal_fails_on_depth(pedigree_path):
    """
    test case for autosomal dominant depth failure
    """

    info_dict = {'gnomad_af': 0.0001, 'gnomad_ac': 0, 'gnomad_hom': 0, 'gene_id': 'TEST1'}
    dom = DominantAutosomal(pedigree=make_flexible_pedigree(pedigree_path))

    # passes with heterozygous
    shallow_variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        depths={'male': 1},
        transcript_consequences=[],
    )
    results = dom.run(principal=shallow_variant)
    assert len(results) == 0


def test_dominant_autosomal_passes(pedigree_path):
    """
    test case for autosomal dominant
    :return:
    """

    info_dict = {'gnomad_af': 0.0001, 'gnomad_ac': 0, 'gnomad_hom': 0, 'cat1': True, 'gene_id': 'TEST1'}

    # attributes relating to categorisation
    boolean_categories = ['cat1']

    dom = DominantAutosomal(pedigree=make_flexible_pedigree(pedigree_path))

    # passes with heterozygous
    passing_variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        depths={'male': 999},
        transcript_consequences=[],
    )
    results = dom.run(principal=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # also passes with homozygous
    passing_variant = SmallVariant(
        info=info_dict,
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        depths={'male': 999},
        transcript_consequences=[],
    )
    results = dom.run(principal=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # no results if no samples
    passing_variant = SmallVariant(
        info=info_dict,
        coordinates=TEST_COORDS,
        boolean_categories=boolean_categories,
        transcript_consequences=[],
    )
    assert len(dom.run(principal=passing_variant)) == 0


@pytest.mark.parametrize('info', [{'gnomad_af': 0.1}, {'gnomad_hom': 2}])
def test_dominant_autosomal_fails(info, pedigree_path):
    """
    test case for autosomal dominant
    :param info: info dict for the variant
    :return:
    """

    dom = DominantAutosomal(pedigree=make_flexible_pedigree(pedigree_path))

    # fails due to high af
    failing_variant = SmallVariant(info=info, het_samples={'male'}, coordinates=TEST_COORDS, transcript_consequences=[])
    assert not dom.run(principal=failing_variant)


def test_recessive_autosomal_hom_passes(pedigree_path):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """
    passing_variant = SmallVariant(
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 1.0},
        depths={'male': 15},
        boolean_categories=['categoryboolean1'],
        info={'categoryboolean1': True, 'gene_id': 'TEST1'},
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=make_flexible_pedigree(pedigree_path))
    results = rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Homozygous'}


def test_recessive_autosomal_hom_passes_with_ab_flag(pedigree_path):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """

    passing_variant = SmallVariant(
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.4},
        depths={'male': 40},
        boolean_categories=['categoryboolean1'],
        info={'categoryboolean1': True, 'gene_id': 'TEST1'},
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=make_flexible_pedigree(pedigree_path))
    results = rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Homozygous'}
    assert passing_variant.get_sample_flags('male') == {'AB Ratio'}


def test_recessive_autosomal_comp_het_male_passes(pedigree_path):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'seqr_link': 'passing1'},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'seqr_link': 'passing2'},
        transcript_consequences=[],
    )
    comp_hets = {'male': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=make_flexible_pedigree(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Comp-Het'}


def test_recessive_autosomal_comp_het_male_passes_partner_flag(pedigree_path):
    """
    info values are defaults (0) & comp-het test is always True we accept a heterozygous variant as a Comp-Het
    """

    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'seqr_link': 'passing1'},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS2,
        ab_ratios={'male': 1.0},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'seqr_link': 'passing2'},
        transcript_consequences=[],
    )
    comp_hets = {'male': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=make_flexible_pedigree(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Comp-Het'}
    assert results[0].flags == {'AB Ratio'}
    assert results[0].support_vars == {'passing2'}


def test_recessive_autosomal_comp_het_female_passes(pedigree_path):
    """
    info values are defaults (0) & comp-het test is always True we accept a heterozygous variant as a Comp-Het
    """
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'seqr_link': 'passing1'},
        transcript_consequences=[],
    )
    passing_variant2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS2,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True, 'seqr_link': 'passing2'},
        transcript_consequences=[],
    )
    comp_hets = {'female': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomalCH(pedigree=make_flexible_pedigree(pedigree_path))
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Comp-Het'}
    assert not results[0].flags
    assert results[0].support_vars == {'passing2'}


def test_recessive_autosomal_comp_het_fails_no_ch_return(pedigree_path):
    """
    check that when the info values are defaults (0) & comp-het test is always False we have no accepted MOI
    """
    failing_variant = SmallVariant(
        info={'gene_id': 'TEST1'},
        het_samples={'male'},
        depths={'male': 50},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalCH(pedigree=make_flexible_pedigree(pedigree_path))
    assert not rec.run(failing_variant)


def test_recessive_autosomal_comp_het_fails_no_paired_call(pedigree_path):
    """
    check that when the info values are defaults (0) & comp-het test is False we have no accepted MOI
    """

    failing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        info={'gene_id': 'TEST1'},
        transcript_consequences=[],
    )
    failing_variant2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS2,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        info={'gene_id': 'TEST1'},
        transcript_consequences=[],
    )

    rec = RecessiveAutosomalCH(pedigree=make_flexible_pedigree(pedigree_path))
    assert not rec.run(
        failing_variant,
        comp_het={'male': {TEST_COORDS2.string_format: [failing_variant2]}},
    )


@pytest.mark.parametrize('info', [{'gnomad_hom': 3, 'gene_id': 'TEST1'}])  # threshold is 2
def test_recessive_autosomal_hom_fails(info, pedigree_path):
    """
    check that when the info values are failures we have no confirmed MOI
    """
    failing_variant = SmallVariant(
        info=info,
        het_samples={'male'},
        hom_samples={'male'},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    rec = RecessiveAutosomalHomo(pedigree=make_flexible_pedigree(pedigree_path))
    assert not rec.run(failing_variant)


def test_x_dominant_female_and_male_het_passes(pedigree_path):
    """
    check that a male and female are accepted as dominant hets
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gnomad_hemi': 0, 'gene_id': 'TEST1', 'categoryboolean1': True},
        het_samples={'female', 'male'},
        depths={'female': 50, 'male': 50},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_dom.run(passing_variant)

    assert len(results) == TWO_EXPECTED
    reasons = {result.reasons.pop() for result in results}
    assert reasons == {'X_Dominant'}


def test_x_dominant_female_hom_passes(pedigree_path):
    """
    check that a female is accepted as a hom
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        hom_samples={'female'},
        depths={'female': 100},
        ab_ratios={'female': 0.5},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant'}


def test_x_dominant_male_hom_passes(pedigree_path):
    """
    check that a male is accepted as a het
    """
    passing_variant = SmallVariant(
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        hom_samples={'male'},
        depths={'male': 100},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XDominant(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant'}


@pytest.mark.parametrize(
    'info,wins',
    [
        ({'gnomad_af': 0.1, 'gene_id': 'TEST1', 'categoryboolean1': True}, 1),
        ({'gnomad_hom': 2, 'gene_id': 'TEST1', 'categoryboolean1': True}, 1),
        ({'gnomad_hemi': 3, 'gene_id': 'TEST1', 'categoryboolean1': False}, 0),
    ],
)
def test_x_dominant_info_fails(info: dict, wins: int, pedigree_path):
    """
    check for info dict exclusions
    """
    passing_variant = SmallVariant(
        info=info,
        hom_samples={'male'},
        coordinates=TEST_COORDS_X_1,
        boolean_categories=['categoryboolean1'],
        transcript_consequences=[],
        depths={'male': 100},
    )
    print(passing_variant)
    x_dom = XDominant(pedigree=make_flexible_pedigree(pedigree_path))
    assert len(x_dom.run(passing_variant)) == wins


def test_x_dominant_female_inactivation_passes(pedigree_path):
    """
    check that a het female, but not a het male, is accepted as dominant
    current test implementation doesn't include family consideration
    """
    info_dict = {
        'gnomad_af': 0.0001,
        'gnomad_ac': 0,
        'gnomad_hom': 0,
        'gene_id': 'TEST1',
        'categoryboolean1': True,
    }
    passing_variant = SmallVariant(
        depths={'female': 100, 'male': 100},
        info=info_dict,
        het_samples={'female', 'male'},
        hom_samples={'male'},
        boolean_categories=['categoryboolean1'],
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    x_dom = XPseudoDominantFemale(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_PseudoDominant'}
    assert 'Affected female with heterozygous variant in XLR gene' in results[0].flags


def test_x_recessive_male_hom_passes(pedigree_path):
    passing_variant = SmallVariant(
        hom_samples={'female', 'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 1.0, 'male': 1.0},
        depths={'female': 100, 'male': 100},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    x_rec = XRecessiveMale(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_rec.run(passing_variant, comp_het={})
    assert len(results) == 1
    assert results[0].reasons == {'X_Male'}


def test_x_recessive_female_hom_passes(pedigree_path):
    """
    :return:
    """

    passing_variant = SmallVariant(
        hom_samples={'female', 'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 1.0, 'male': 1.0},
        depths={'female': 100, 'male': 100},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    x_rec = XRecessiveFemaleHom(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_rec.run(passing_variant, comp_het={})
    assert len(results) == 1
    assert results[0].reasons == {'X_Recessive HOM Female'}


def test_x_recessive_male_het_passes(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        boolean_categories=['categoryboolean1'],
        info={'gene_id': 'TEST1', 'categoryboolean1': True},
        transcript_consequences=[],
    )
    x_rec = XRecessiveMale(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Male'}


def test_x_recessive_female_het_passes(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'seqr_link': 'passing1'},
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['female'], 'seqr_link': 'passing2'},
        transcript_consequences=[],
    )
    comp_hets = {'female': {'X-1-G-T': [passing_variant_2]}}
    x_rec = XRecessiveFemaleCH(pedigree=make_flexible_pedigree(pedigree_path))
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'X_RecessiveFemaleCompHet'}
    assert results[0].support_vars == {'passing2'}


def test_het_de_novo_passes(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        sample_categories=['categorysample4'],
        ab_ratios={'female': 0.5},
        depths={'female': 99},
        info={'gene_id': 'TEST1', 'categorysample4': ['female']},
        transcript_consequences=[],
    )
    dom_a = DominantAutosomal(pedigree=make_flexible_pedigree(pedigree_path))
    results = dom_a.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}
    assert not results[0].flags


def test_het_de_novo_het_passes_flagged(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        sample_categories=['categorysample4'],
        ab_ratios={'female': 0.5},
        depths={'female': 99},
        info={'gene_id': 'TEST1', 'categorysample4': ['female']},
        transcript_consequences=[],
    )
    dom_a = DominantAutosomal(pedigree=make_flexible_pedigree(pedigree_path))
    results = dom_a.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}


def test_x_recessive_female_het_fails(pedigree_path):
    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['male']},
        transcript_consequences=[],
    )
    passing_variant_2 = SmallVariant(
        het_samples={'male'},
        coordinates=TEST_COORDS_X_2,
        ab_ratios={'male': 0.5},
        depths={'male': 50},
        sample_categories=['categorysample4'],
        info={'gene_id': 'TEST1', 'categorysample4': ['male']},
        transcript_consequences=[],
    )
    comp_hets = {'female': {'x-2-A-C': [passing_variant_2]}}
    x_rec = XRecessiveFemaleCH(pedigree=make_flexible_pedigree(pedigree_path))
    assert not x_rec.run(passing_variant, comp_het=comp_hets)


@mock.patch('talos.moi_tests.check_for_second_hit')
def test_x_recessive_female_het_no_pair_fails(second_hit: mock.Mock, pedigree_path):
    """ """

    passing_variant = SmallVariant(
        het_samples={'female'},
        coordinates=TEST_COORDS_X_1,
        ab_ratios={'female': 0.5},
        depths={'female': 50},
        info={'gene_id': 'TEST1', 'categorysample1': True, 'boolean_categories': 'categorysample1'},
        transcript_consequences=[],
    )
    second_hit.return_value = []
    assert not XRecessiveFemaleCH(pedigree=make_flexible_pedigree(pedigree_path)).run(passing_variant)


def test_check_familial_inheritance_simple(pedigree_path):
    """
    test the check_familial_inheritance method
    trio male, mother_1, father_1; only 'male' is affected
    """

    base_moi = BaseMoi(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')
    assert base_moi.check_familial_inheritance(sample_id='male', called_variants={'male'})


def test_check_familial_inheritance_mother_fail(pedigree_path):
    """
    test the check_familial_inheritance method
    """

    base_moi = BaseMoi(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')
    assert not base_moi.check_familial_inheritance(sample_id='male', called_variants={'male', 'mother_1'})


def test_check_familial_inheritance_mother_passes(pedigree_path):
    """
    test the check_familial_inheritance method
    mother in variant calls, but partial penetrance
    """

    base_moi = BaseMoi(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')

    assert base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'mother_1'},
        partial_pen=True,
    )


def test_check_familial_inheritance_father_fail(pedigree_path):
    """
    test the check_familial_inheritance method
    """

    base_moi = BaseMoi(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')
    assert not base_moi.check_familial_inheritance(sample_id='male', called_variants={'male', 'father_1'})


def test_check_familial_inheritance_father_passes(pedigree_path):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    """

    base_moi = BaseMoi(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'father_1'},
        partial_pen=True,
    )
    assert result


def test_check_familial_inheritance_top_down(pedigree_path):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    """

    base_moi = BaseMoi(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')
    assert base_moi.check_familial_inheritance(
        sample_id='father_1',
        called_variants={'male', 'father_1'},
        partial_pen=True,
    )


def test_check_familial_inheritance_no_calls(pedigree_path):
    """
    test the check_familial_inheritance method where there are no calls
    we lazily pass this - the assumption is that we're only assessing samples with variants
    this method checks parents as a trio, not the sample in question
    """

    base_moi = BaseMoi(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')
    assert base_moi.check_familial_inheritance(sample_id='male', called_variants=set(), partial_pen=True)


def test_genotype_calls(pedigree_path):
    """
    test the manual genotype assignments
    """
    base_moi = DominantAutosomal(pedigree=make_flexible_pedigree(pedigree_path), applied_moi='applied')

    info_dict = {'gnomad_af': 0.0001, 'gnomad_ac': 0, 'gnomad_hom': 0, 'gene_id': 'TEST1'}
    variant = SmallVariant(
        info=info_dict,
        het_samples={'male'},
        hom_samples={'female'},
        coordinates=TEST_COORDS,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(variant, 'male') == {'father_1': 'WT', 'male': 'Het', 'mother_1': 'WT'}
    assert base_moi.get_family_genotypes(variant, 'female') == {'father_2': 'WT', 'female': 'Hom', 'mother_2': 'WT'}
    x_variant = SmallVariant(
        info=info_dict,
        het_samples={'male', 'female'},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(x_variant, 'male') == {'father_1': 'WT', 'male': 'Hemi', 'mother_1': 'WT'}
    assert base_moi.get_family_genotypes(x_variant, 'female') == {'father_2': 'WT', 'female': 'Het', 'mother_2': 'WT'}

    x_variant_2 = SmallVariant(
        info=info_dict,
        hom_samples={'male', 'female'},
        coordinates=TEST_COORDS_X_1,
        transcript_consequences=[],
    )
    assert base_moi.get_family_genotypes(x_variant_2, 'male') == {'father_1': 'WT', 'male': 'Hemi', 'mother_1': 'WT'}
    assert base_moi.get_family_genotypes(x_variant_2, 'female') == {'father_2': 'WT', 'female': 'Hom', 'mother_2': 'WT'}
    variant_missing = SmallVariant(info=info_dict, coordinates=TEST_COORDS, transcript_consequences=[])
    assert base_moi.get_family_genotypes(variant_missing, 'male') == {'father_1': 'WT', 'male': 'WT', 'mother_1': 'WT'}
