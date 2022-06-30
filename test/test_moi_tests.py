"""
tests relating to the MOI filters
"""


from dataclasses import dataclass, field
from typing import Any, Dict, List, Set

from unittest import mock

import pytest
from reanalysis.moi_tests import (
    check_for_second_hit,
    BaseMoi,
    DominantAutosomal,
    GNOMAD_AD_AC_THRESHOLD,
    GNOMAD_DOM_HOM_THRESHOLD,
    GNOMAD_HEMI_THRESHOLD,
    GNOMAD_RARE_THRESHOLD,
    GNOMAD_REC_HOM_THRESHOLD,
    MOIRunner,
    RecessiveAutosomal,
    XDominant,
    XRecessive,
)

from reanalysis.utils import Coordinates


MOI_CONF = {
    GNOMAD_REC_HOM_THRESHOLD: 2,
    GNOMAD_DOM_HOM_THRESHOLD: 1,
    GNOMAD_AD_AC_THRESHOLD: 10,
    GNOMAD_RARE_THRESHOLD: 0.01,
    GNOMAD_HEMI_THRESHOLD: 2,
}
TEST_COORDS = Coordinates('1', 1, 'A', 'C')
TEST_COORDS2 = Coordinates('2', 2, 'G', 'T')
TINY_CONFIG = {'male': 'male'}
TINY_COMP_HET = {}


@dataclass
class SimpleVariant:
    """
    a fake version of AbstractVariant
    """

    info: Dict[str, Any]
    het_samples: Set[str]
    hom_samples: Set[str]
    coords: Coordinates
    category_1: bool = True
    category_4: List[str] = field(default_factory=list)

    def sample_specific_category_check(self, sample):
        """
        pass
        :param sample:
        :return:
        """
        return self.category_1 or sample in self.category_4


@dataclass
class RecessiveSimpleVariant:
    """
    a fake version of AbstractVariant
    """

    info: Dict[str, Any]
    het_samples: Set[str]
    hom_samples: Set[str]
    coords: Coordinates
    category_4: List[str]
    # add category default
    category_1: bool = True

    @property
    def category_1_2_3(self):
        """
        mock method
        :return:
        """
        return self.category_1

    def sample_de_novo(self, sample):
        """
        pass
        :param sample:
        :return:
        """
        return sample in self.category_4

    def sample_specific_category_check(self, sample):
        """
        pass
        :param sample:
        :return:
        """
        return (sample in self.category_4) or self.category_1_2_3


@pytest.mark.parametrize(
    'first,comp_hets,sample,values',
    (
        ('', {}, '', []),  # no values
        ('', {}, 'a', []),  # sample not present
        ('', {'a': {'foo': []}}, 'a', []),  # var not present
        (
            'foo',
            {'a': {'foo': ['bar']}},
            'a',
            ['bar'],
        ),  # all values present
        (
            'foo',
            {'a': {'foo': ['bar', 'baz']}},
            'a',
            ['bar', 'baz'],
        ),  # all values present
    ),
)
def test_check_second_hit(first, comp_hets, sample, values):
    """
    quick test for the 2nd hit mechanic
    return all strings when the comp-het lookup contains:
        - the sample
        - the gene
        - the variant signature
    :return:
    """

    assert (
        check_for_second_hit(first_variant=first, comp_hets=comp_hets, sample=sample)
        == values
    )


@pytest.mark.parametrize(
    'moi_string,filters',
    (
        ('Monoallelic', ['DominantAutosomal']),
        ('Mono_And_Biallelic', ['DominantAutosomal', 'RecessiveAutosomal']),
        ('Unknown', ['DominantAutosomal', 'RecessiveAutosomal']),
        ('Biallelic', ['RecessiveAutosomal']),
        (
            'Hemi_Mono_In_Female',
            ['XRecessive', 'XDominant'],
        ),
        ('Hemi_Bi_In_Female', ['XRecessive']),
        ('Y_Chrom_Variant', ['YHemi']),
    ),
)
def test_moi_runner(moi_string: str, filters: List[str], peddy_ped):
    """

    :param moi_string:
    :param filters:
    :return:
    """
    test_runner = MOIRunner(
        pedigree=peddy_ped,
        target_moi=moi_string,
        config=TINY_CONFIG,
    )

    # string-comparison
    # the imported (uninstantiated) objects don't have __class__
    # and the instantiated objects don't have a __name__
    for filter1, filter2 in zip(test_runner.filter_list, filters):
        assert filter2 in str(filter1.__class__)


def test_dominant_autosomal_passes(peddy_ped):
    """
    test case for autosomal dominant
    :return:
    """

    info_dict = {'gnomad_af': 0.0001, 'gnomad_ac': 0, 'gnomad_hom': 0}

    dom = DominantAutosomal(pedigree=peddy_ped, config=MOI_CONF)

    # passes with heterozygous
    passing_variant = SimpleVariant(
        info=info_dict, het_samples={'male'}, hom_samples=set(), coords=TEST_COORDS
    )
    results = dom.run(principal_var=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # also passes with homozygous
    passing_variant = SimpleVariant(
        info=info_dict, het_samples=set(), hom_samples={'male'}, coords=TEST_COORDS
    )
    results = dom.run(principal_var=passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Dominant'}

    # no results if no samples
    passing_variant = SimpleVariant(
        info=info_dict, het_samples=set(), hom_samples=set(), coords=TEST_COORDS
    )
    assert not dom.run(principal_var=passing_variant)


@pytest.mark.parametrize(
    'info',
    [{'gnomad_af': 0.1}, {'gnomad_hom': 2}],
)
def test_dominant_autosomal_fails(info, peddy_ped):
    """
    test case for autosomal dominant
    :param info: info dict for the variant
    :return:
    """

    dom = DominantAutosomal(pedigree=peddy_ped, config=MOI_CONF)

    # fails due to high af
    failing_variant = SimpleVariant(
        info=info, het_samples={'male'}, hom_samples=set(), coords=TEST_COORDS
    )
    assert not dom.run(principal_var=failing_variant)


def test_recessive_autosomal_hom_passes(peddy_ped):
    """
    check that when the info values are defaults (0)
    we accept a homozygous variant as a Recessive
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples=set(),
        hom_samples={'male'},
        coords=TEST_COORDS,
        category_4=[],
    )
    rec = RecessiveAutosomal(pedigree=peddy_ped, config={GNOMAD_REC_HOM_THRESHOLD: 1})
    results = rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Homozygous'}


def test_recessive_autosomal_comp_het_male_passes(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=TEST_COORDS,
        category_4=[],
    )
    passing_variant2 = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=TEST_COORDS2,
        category_4=[],
    )
    comp_hets = {'male': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomal(pedigree=peddy_ped, config={GNOMAD_REC_HOM_THRESHOLD: 1})
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Compound-Het'}


def test_recessive_autosomal_comp_het_female_passes(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always True
    we accept a heterozygous variant as a Comp-Het
    :return:
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=TEST_COORDS,
        category_4=[],
    )
    passing_variant2 = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=TEST_COORDS2,
        category_4=[],
    )
    comp_hets = {'female': {TEST_COORDS.string_format: [passing_variant2]}}
    rec = RecessiveAutosomal(pedigree=peddy_ped, config={GNOMAD_REC_HOM_THRESHOLD: 1})
    results = rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'Autosomal Recessive Compound-Het'}


def test_recessive_autosomal_comp_het_fails_no_ch_return(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always False
    we have no accepted MOI

    :return:
    """

    failing_variant = SimpleVariant(
        info={}, het_samples={'male'}, hom_samples=set(), coords=TEST_COORDS
    )
    rec = RecessiveAutosomal(pedigree=peddy_ped, config={GNOMAD_REC_HOM_THRESHOLD: 1})
    assert not rec.run(failing_variant)


def test_recessive_autosomal_comp_het_fails_no_paired_call(peddy_ped):
    """
    check that when the info values are defaults (0)
    and the comp-het test is always False
    we have no accepted MOI

    :return:
    """

    failing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=TEST_COORDS,
        category_4=[],
    )
    failing_variant2 = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=TEST_COORDS2,
        category_4=[],
    )

    rec = RecessiveAutosomal(pedigree=peddy_ped, config={GNOMAD_REC_HOM_THRESHOLD: 1})
    assert not rec.run(
        failing_variant,
        comp_het={'male': {TEST_COORDS2.string_format: [failing_variant2]}},
    )


@pytest.mark.parametrize(
    'info',
    [{'gnomad_hom': 2}],
)
def test_recessive_autosomal_hom_fails(info, peddy_ped):
    """
    check that when the info values are failures
    we have no confirmed MOI
    """

    failing_variant = SimpleVariant(
        info=info, het_samples={'male'}, hom_samples={'male'}, coords=TEST_COORDS
    )
    rec = RecessiveAutosomal(pedigree=peddy_ped, config={GNOMAD_REC_HOM_THRESHOLD: 1})
    assert not rec.run(failing_variant)


def test_x_dominant_female_and_male_het_passes(peddy_ped):
    """
    check that a male is accepted as a het
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info={'gnomad_hemi': 0},
        het_samples={'female', 'male'},
        hom_samples=set(),
        coords=x_coords,
    )
    x_dom = XDominant(pedigree=peddy_ped, config=MOI_CONF)
    results = x_dom.run(passing_variant)

    assert len(results) == 2
    reasons = sorted([result.reasons.pop() for result in results])
    assert reasons == ['X_Dominant Female', 'X_Dominant Male']


def test_x_dominant_female_hom_passes(peddy_ped):
    """
    check that a male is accepted as a het
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info={'gnomad_hemi': 0},
        hom_samples={'female'},
        het_samples=set(),
        coords=x_coords,
    )
    x_dom = XDominant(pedigree=peddy_ped, config=MOI_CONF)
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant Female'}


def test_x_dominant_male_hom_passes(peddy_ped):
    """
    check that a male is accepted as a het
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info={'gnomad_hemi': 0},
        hom_samples={'male'},
        het_samples=set(),
        coords=x_coords,
    )
    x_dom = XDominant(pedigree=peddy_ped, config=MOI_CONF)
    results = x_dom.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Dominant Male'}


@pytest.mark.parametrize(
    'info',
    [
        {'gnomad_af': 0.1},
        {'gnomad_hom': 2},
        {'gnomad_hemi': 3},
    ],
)
def test_x_dominant_info_fails(info, peddy_ped):
    """
    check for info dict exclusions
    :param info:
    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = SimpleVariant(
        info=info, hom_samples={'male'}, het_samples=set(), coords=x_coords
    )
    x_dom = XDominant(pedigree=peddy_ped, config=MOI_CONF)
    assert not x_dom.run(passing_variant)


def test_x_recessive_male_and_female_hom_passes(peddy_ped):
    """

    :return:
    """

    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = RecessiveSimpleVariant(
        info={},
        hom_samples={'female', 'male'},
        het_samples=set(),
        coords=x_coords,
        category_4=[],
    )
    x_rec = XRecessive(pedigree=peddy_ped, config=MOI_CONF)
    results = x_rec.run(passing_variant, comp_het={})
    assert len(results) == 2

    reasons = sorted([result.reasons.pop() for result in results])
    assert reasons == ['X_Recessive Female', 'X_Recessive Male']


def test_x_recessive_male_het_passes(peddy_ped):
    """

    :return:
    """
    x_coords = Coordinates('x', 1, 'A', 'C')
    passing_variant = RecessiveSimpleVariant(
        info={}, het_samples={'male'}, hom_samples=set(), coords=x_coords, category_4=[]
    )
    x_rec = XRecessive(pedigree=peddy_ped, config=MOI_CONF)
    results = x_rec.run(passing_variant)
    assert len(results) == 1
    assert results[0].reasons == {'X_Recessive Male'}


def test_x_recessive_female_het_passes(peddy_ped):
    """

    :return:
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 1, 'A', 'C'),
        category_4=['female'],
    )
    passing_variant_2 = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 2, 'A', 'C'),
        category_4=['female'],
    )
    comp_hets = {'female': {'x-1-A-C': [passing_variant_2]}}
    x_rec = XRecessive(pedigree=peddy_ped, config=MOI_CONF)
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    assert len(results) == 1
    assert results[0].reasons == {'X_Recessive Compound-Het Female'}


def test_x_recessive_female_het_fails(peddy_ped):
    """

    :return:
    """

    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 1, 'A', 'C'),
        category_4=['male'],
    )
    passing_variant_2 = RecessiveSimpleVariant(
        info={},
        het_samples={'male'},
        hom_samples=set(),
        coords=Coordinates('x', 2, 'A', 'C'),
        category_4=['male'],
    )
    comp_hets = {'female': {'x-2-A-C': [passing_variant_2]}}
    x_rec = XRecessive(pedigree=peddy_ped, config=MOI_CONF)
    results = x_rec.run(passing_variant, comp_het=comp_hets)
    assert not results


@mock.patch('reanalysis.moi_tests.check_for_second_hit')
def test_x_recessive_female_het_no_pair_fails(second_hit: mock.patch, peddy_ped):
    """

    :return:
    """

    second_hit.return_value = []
    passing_variant = RecessiveSimpleVariant(
        info={},
        het_samples={'female'},
        hom_samples=set(),
        coords=Coordinates('x', 1, 'A', 'C'),
        category_4=[],
    )
    x_rec = XRecessive(pedigree=peddy_ped, config=MOI_CONF)
    assert not x_rec.run(passing_variant)


# trio male, mother_1, father_1; only 'male' is affected
def test_check_familial_inheritance_simple(peddy_ped):
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, config=TINY_CONFIG, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male'}
    )
    assert result


def test_check_familial_inheritance_mother_fail(peddy_ped):
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, config=TINY_CONFIG, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male', 'mother_1'}
    )
    assert not result


def test_check_familial_inheritance_mother_passes(peddy_ped):
    """
    test the check_familial_inheritance method
    mother in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, config=TINY_CONFIG, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'mother_1'},
        partial_penetrance=True,
    )
    assert result


def test_check_familial_inheritance_father_fail(peddy_ped):
    """
    test the check_familial_inheritance method
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, config=TINY_CONFIG, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male', called_variants={'male', 'father_1'}
    )
    assert not result


def test_check_familial_inheritance_father_passes(peddy_ped):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, config=TINY_CONFIG, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants={'male', 'father_1'},
        partial_penetrance=True,
    )
    assert result


def test_check_familial_inheritance_top_down(peddy_ped):
    """
    test the check_familial_inheritance method
    father in variant calls, but partial penetrance
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, config=TINY_CONFIG, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='father_1',
        called_variants={'male', 'father_1'},
        partial_penetrance=True,
    )
    assert result


def test_check_familial_inheritance_no_calls(peddy_ped):
    """
    test the check_familial_inheritance method where there are no calls
    will fail as affected proband not in calls
    :return:
    """

    base_moi = BaseMoi(pedigree=peddy_ped, config=TINY_CONFIG, applied_moi='applied')

    result = base_moi.check_familial_inheritance(
        sample_id='male',
        called_variants=set(),
        partial_penetrance=True,
    )
    # should fail immediately
    assert not result
