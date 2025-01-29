"""
A number of classes, each representing one Mode of Inheritance
One class (MoiRunner) to run all the appropriate MOIs on a variant

Reduce the PanelApp plain text MOI description into a few categories
We then run a permissive MOI match for the variant
"""

# mypy: ignore-errors
from abc import abstractmethod

from talos.config import config_retrieve
from talos.models import VARIANT_MODELS, Pedigree, ReportVariant, SmallVariant, StructuralVariant
from talos.utils import X_CHROMOSOME, CompHetDict

# config keys to use for dominant MOI tests
CALLSET_AF_SV_DOMINANT = 'callset_af_sv_dominant'
CB1 = 'categoryboolean1'
GNOMAD_RARE_THRESHOLD = 'gnomad_dominant'
GNOMAD_AD_AC_THRESHOLD = 'gnomad_max_ac_dominant'
GNOMAD_DOM_HOM_THRESHOLD = 'gnomad_max_homs_dominant'
GNOMAD_REC_HOM_THRESHOLD = 'gnomad_max_homs_recessive'
GNOMAD_HEMI_THRESHOLD = 'gnomad_max_hemi'
INFO_HOMS = {'gnomad_hom', 'gnomad_ex_hom'}
INFO_HEMI = {'gnomad_hemi', 'gnomad_ex_hemi'}
SV_AF_KEY = 'gnomad_v2.1_sv_AF'
SV_HEMI = {'male_n_hemialt'}
SV_HOMS = {'male_n_homalt', 'female_n_homalt'}


def check_for_second_hit(
    first_variant: str,
    comp_hets: CompHetDict,
    sample: str,
    require_non_support: bool = False,
) -> list[VARIANT_MODELS]:
    """
    checks for a second hit partner in this gene

    Example formatting of the comp-het dict
    {
        "SampleID": {
            "12-52287177-T-C": [
                Variant(12-52287180-TGG-T)
            ],
            "12-52287180-TGG-T": [
                Variant(12-52287177-T-C)
            ]
        } ...
    }

    Args:
        first_variant (str): string representation of variant1
        comp_hets (dict[str, Variant]): lookup for compound hets
        sample (str): sample ID
        require_non_support (bool): if true, don't return Support only

    Returns:
        a list of variants which are potential partners
    """

    # check if the sample has any comp-hets
    if sample not in comp_hets:
        return []

    partners = comp_hets[sample].get(first_variant, [])
    if require_non_support:
        return [partner for partner in partners if not partner.sample_support_only(sample)]
    return partners


class MOIRunner:
    """
    The abstract class for a single MOI runner
    This will be instantiated once per MOI, and run once per related gene, on the collection of all variants in the gene
    """

    def __init__(self, pedigree: Pedigree, target_moi: str):
        """
        for each possible MOI, choose the appropriate filters to apply
        ran into a situation where the ID of target_moi didn't match the
        exact same MOI as the IDs were different.

        This logic is only called once per MOI, not once per variant

        Args:
            pedigree ():
            target_moi ():
        """

        # for unknown, we catch all possible options?
        # should we be doing both checks for Monoallelic?
        if target_moi == 'Monoallelic':
            self.filter_list = [DominantAutosomal(pedigree=pedigree)]
        elif target_moi in ['Mono_And_Biallelic', 'Unknown']:
            self.filter_list = [
                DominantAutosomal(pedigree=pedigree),
                RecessiveAutosomalHomo(pedigree=pedigree),
                RecessiveAutosomalCH(pedigree=pedigree),
            ]
        elif target_moi == 'Biallelic':
            self.filter_list = [RecessiveAutosomalHomo(pedigree=pedigree), RecessiveAutosomalCH(pedigree=pedigree)]

        elif target_moi == 'Hemi_Mono_In_Female':
            self.filter_list = [XRecessiveMale(pedigree=pedigree), XDominant(pedigree=pedigree)]

        elif target_moi == 'Hemi_Bi_In_Female':
            self.filter_list = [
                XRecessiveMale(pedigree=pedigree),
                XRecessiveFemaleHom(pedigree=pedigree),
                XRecessiveFemaleCH(pedigree=pedigree),
                XPseudoDominantFemale(pedigree=pedigree),
            ]

        else:
            raise KeyError(f'MOI type {target_moi} is not addressed in MOI')

    def run(self, principal_var, comp_het: CompHetDict | None = None, partial_pen: bool = False) -> list[ReportVariant]:
        """
        run method - triggers each relevant inheritance model

        Args:
            principal_var (): the variant we are focused on
            comp_het ():
            partial_pen ():
        """

        if comp_het is None:
            comp_het = {}

        moi_matched = []
        for model in self.filter_list:
            moi_matched.extend(model.run(principal=principal_var, comp_het=comp_het, partial_pen=partial_pen))
        return moi_matched


class BaseMoi:
    """
    Definition of the MOI base class
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str):
        """
        base class
        """
        if applied_moi is None:
            raise ValueError('An applied MOI needs to reach the Base Class')
        self.pedigree = pedigree
        self.applied_moi = applied_moi
        self.minimum_depth = config_retrieve(['RunHailFiltering', 'minimum_depth'], 10)

    @abstractmethod
    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        run all applicable inheritance patterns and finds good fits
        """

    def check_affected_category_depths(self, variant: VARIANT_MODELS, sample_id: str) -> bool:
        """

        Args:
            variant ():
            sample_id ():

        Returns:
            bool: True if the variant passes the depth checks
        """
        if (
            not (
                self.pedigree.by_id[sample_id].affected == '2'
                and variant.sample_category_check(sample_id, allow_support=False)
            )
        ) or variant.check_read_depth(sample_id, self.minimum_depth, var_is_cat_1=variant.info.get(CB1)):
            return True
        return False

    def check_familial_inheritance(self, sample_id: str, called_variants: set[str], partial_pen: bool = False) -> bool:
        """
        check for single variant inheritance
        - find the family ID from this sample
        - iterate through all family members, and check that MOI holds for all

        Return False if a participant fails tests, True if all checked are ok

        NOTE: this called_variants pool is prepared before calling this method.
        If we only want to check hom calls, only send hom calls. If we are checking for
        dominant conditions, we would bundle hom and het calls into this set

        Sex chrom MOI is being broken down in a more granular way, so this may need
        revisiting

        Args:
            sample_id ():
            called_variants (set[str]): the set of sample_ids which have this variant
            partial_pen (bool): if True, permit unaffected has variant call

        Returns:
            True if all tests pass, else False
        """

        this_member = self.pedigree.by_id[sample_id]

        # check for valid inheritance within the immediate trio, if possible
        for member_id in [this_member.father, this_member.mother]:
            if member_id is None or (iter_member := self.pedigree.by_id.get(member_id)) is None:
                continue

            # complete & incomplete penetrance - affected samples must have the variant
            # complete pen. requires participants to be affected if they have the var
            # if any of these combinations occur, fail the family
            if (iter_member.affected == '2' and member_id not in called_variants) or (
                member_id in called_variants and not partial_pen and iter_member.affected != '2'
            ):
                return False

        return True

    def get_family_genotypes(self, variant: VARIANT_MODELS, sample_id: str) -> dict[str, str]:
        """

        Args:
            variant (VARIANT_MODELS):
            sample_id (str): the sample ID to gather genotypes for

        Returns:
            list[str]: a list of all the participants, and GTs
        """

        def get_sample_genotype(member_id: str, sex: str) -> str:
            """
            for this specific member, find the genotype
            Args:
                member_id (str): sample ID in the pedigree
                sex (str): male/female/unknown

            Returns:
                str: text representation of this genotype
            """

            if variant.coordinates.chrom in X_CHROMOSOME:
                if sex == '1' and (member_id in variant.het_samples or member_id in variant.hom_samples):
                    return 'Hemi'

                if member_id in variant.het_samples:
                    return 'Het'
                if member_id in variant.hom_samples:
                    return 'Hom'

            elif member_id in variant.het_samples:
                return 'Het'
            elif member_id in variant.hom_samples:
                return 'Hom'

            return 'WT'

        sample_family_id = self.pedigree.by_id[sample_id].family
        return {
            member.id: get_sample_genotype(member_id=member.id, sex=member.sex)
            for member in self.pedigree.by_family[sample_family_id]
        }

    @staticmethod
    def check_frequency_fails(info: dict, thresholds: dict[str, int | float], permit_clinvar: bool = True) -> bool:
        """
        Method to check multiple info keys against a single threshold
        This just reduces the line count, as this is called a bunch of times
        By default we always let clinvar pathogenic variants pass

        Args:
            info (): the dict of values for this dict
            thresholds (): the dict of keys - thresholds to test against
            permit_clinvar (bool): if True, always allow clinvar pathogenic variants to pass

        Returns:
            True if any of the info attributes are above the threshold, unless we use a clinvar escape
        """
        if permit_clinvar and info.get(CB1):
            return False
        return any(info.get(key, 0) > test for key, test in thresholds.items())

    @staticmethod
    def check_callset_af_fails(info: dict, threshold: float) -> bool:
        """
        unlike the frequency method above, we want to apply an escape condition here
        - if the callset is large enough, apply this filter

        this is predicated on info containing both ac and af
        previously we've permitted < 5 instances in the callset through

        Args:
            info (dict): info dict for this variant
            threshold (float): the threshold value for this test

        Returns:
            True if this variant is above the filtering threshold
        """

        min_ac = 5
        if info.get('ac', 0) <= min_ac:
            return False

        return info.get('af', 0.0) >= threshold

    def check_comp_het(self, sample_id: str, variant_1: VARIANT_MODELS, variant_2: VARIANT_MODELS) -> bool:
        """
        use parents to accept or dismiss the comp-het
        If the 'comp-het' pair are inherited from a single parent, they are in cis
        rather than trans, and reporting as a comp-het would be misleading.

        compound het is inherently not inherited from a single parent, so rule out
        when either parent has both, or either parent is affected

        Args:
            sample_id (str): sample ID to check for
            variant_1 (VARIANT_MODELS): first variant of comp-het pair
            variant_2 (VARIANT_MODELS): second variant of comp-het pair

        Returns:
            bool: True if these two variants form a comp-het
        """

        # if both vars are present in a single parent: not a compound het
        # or if the parent is affected: not causative
        sample_ped_entry = self.pedigree.by_id[sample_id]
        for parent in [sample_ped_entry.mother, sample_ped_entry.father]:
            # skip to prevent crashing on !trios
            if parent is None:
                continue

            if ((parent in variant_1.het_samples) and (parent in variant_2.het_samples)) or self.pedigree.by_id[
                parent
            ].affected == '2':
                return False
        return True


class DominantAutosomal(BaseMoi):
    """
    This class can also be called by the X-linked Dominant, in which case the
    Applied_MOI by name is overridden
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'Autosomal Dominant'):
        """
        Simplest: AD MOI
        """

        self.ad_threshold = config_retrieve(['ValidateMOI', GNOMAD_RARE_THRESHOLD])
        self.ac_threshold = config_retrieve(['ValidateMOI', GNOMAD_AD_AC_THRESHOLD])
        self.hom_threshold = config_retrieve(['ValidateMOI', GNOMAD_DOM_HOM_THRESHOLD])
        self.sv_af_threshold = config_retrieve(['ValidateMOI', CALLSET_AF_SV_DOMINANT])

        # prepare the AF test dicts
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS}
            | {'gnomad_ac': self.ac_threshold, 'gnomad_af': self.ad_threshold},
            StructuralVariant.__name__: {'af': self.sv_af_threshold, SV_AF_KEY: self.sv_af_threshold},
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        Simplest MOI, exclusions based on HOM count and AF
        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # reject support for dominant MOI, apply checks based on var type
        if (
            principal.support_only
            or self.check_frequency_fails(
                principal.info,
                self.freq_tests[principal.__class__.__name__],
            )
            or self.check_callset_af_fails(principal.info, self.ad_threshold)
        ):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        samples_with_this_variant = principal.het_samples.union(principal.hom_samples)
        for sample_id in samples_with_this_variant:
            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised
            # force a minimum depth on the proband call
            if not (
                self.pedigree.by_id[sample_id].affected == '2'
                and principal.sample_category_check(sample_id, allow_support=False)
            ) or (
                principal.check_read_depth(
                    sample_id,
                    self.minimum_depth,
                    var_is_cat_1=principal.info.get(CB1),
                )
            ):
                continue

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree.by_id[sample_id].family,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
                ),
            )

        return classifications


class RecessiveAutosomalCH(BaseMoi):
    """
    inheritance test class for Recessive inheritance
    requires single hom variant, or compound het
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'Autosomal Recessive Comp-Het'):
        """ """
        self.hom_threshold = config_retrieve(['ValidateMOI', GNOMAD_REC_HOM_THRESHOLD])
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {key: self.hom_threshold for key in SV_HOMS},
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,  # noqa: ARG002
    ) -> list[ReportVariant]:
        """
        valid if present as compound het
        counts as being phased if a compound het is split between parents

        Args:
            principal (VARIANT_MODELS): main variant being evaluated
            comp_het (dict): comp-het partners
            partial_pen (bool):

        Returns:
            list[ReportVariant]: data object if RecessiveAutosomal fits
        """

        if comp_het is None:
            comp_het = {}

        classifications = []

        # remove if too many homs are present in population databases
        if self.check_frequency_fails(
            principal.info,
            self.freq_tests[principal.__class__.__name__],
        ):
            return classifications

        # if hets are present, try and find support
        for sample_id in principal.het_samples:
            # skip primary analysis for unaffected members
            # this sample must be categorised - check Cat 4 contents
            if (
                not (
                    self.pedigree.by_id[sample_id].affected == '2'
                    and principal.sample_category_check(sample_id, allow_support=True)
                )
            ) or (principal.check_read_depth(sample_id, self.minimum_depth, principal.info.get(CB1))):
                continue

            for partner_variant in check_for_second_hit(
                first_variant=principal.coordinates.string_format,
                comp_hets=comp_het,
                sample=sample_id,
                require_non_support=principal.sample_support_only(sample_id),
            ):
                if partner_variant.check_read_depth(
                    sample_id,
                    self.minimum_depth,
                    partner_variant.info.get(CB1),
                ) or self.check_frequency_fails(
                    partner_variant.info,
                    self.freq_tests[partner_variant.__class__.__name__],
                ):
                    continue

                # check if this is a candidate for comp-het inheritance
                if not self.check_comp_het(sample_id=sample_id, variant_1=principal, variant_2=partner_variant):
                    continue

                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.by_id[sample_id].family,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        reasons={self.applied_moi},
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        support_vars={partner_variant.info['seqr_link']},
                        flags=principal.get_sample_flags(sample_id) | partner_variant.get_sample_flags(sample_id),
                        independent=False,
                    ),
                )

        return classifications


class RecessiveAutosomalHomo(BaseMoi):
    """
    inheritance test class for Recessive inheritance
    requires single hom variant, or compound het
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'Autosomal Recessive Homozygous'):
        """ """
        self.hom_threshold = config_retrieve(['ValidateMOI', GNOMAD_REC_HOM_THRESHOLD])
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {key: self.hom_threshold for key in SV_HOMS},
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        explicitly tests HOMs

        Args:
            principal (VARIANT_MODELS): main variant being evaluated
            comp_het (dict): comp-het partners
            partial_pen (bool):

        Returns:
            list[ReportVariant]: data object if RecessiveAutosomal fits
        """

        classifications = []

        # remove if too many homs are present in population databases
        if principal.support_only or self.check_frequency_fails(
            principal.info,
            self.freq_tests[principal.__class__.__name__],
        ):
            return classifications

        for sample_id in principal.hom_samples:
            # skip primary analysis for unaffected members
            # require this sample to be categorised - check Sample contents
            # minimum depth of call
            if (
                not (
                    self.pedigree.by_id[sample_id].affected == '2'
                    and principal.sample_category_check(sample_id, allow_support=False)
                )
            ) or principal.check_read_depth(sample_id, self.minimum_depth, principal.info.get(CB1)):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=principal.hom_samples,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree.by_id[sample_id].family,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
                ),
            )

        return classifications


class XDominant(BaseMoi):
    """
    for males and females, accept het
    effectively the same as DominantAutosomal?
    just override the type, and use AD

    GATK MALES ARE CALLED HOM (OR HET pseudo-autosomal?)
    re-implement here, but don't permit Male X-Homs
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'X_Dominant'):
        """
        accept male hets and homs, and female hets without support

        Args:
            pedigree ():
            applied_moi ():
        """
        self.ad_threshold = config_retrieve(['ValidateMOI', GNOMAD_RARE_THRESHOLD])
        self.ac_threshold = config_retrieve(['ValidateMOI', GNOMAD_AD_AC_THRESHOLD])
        self.hom_threshold = config_retrieve(['ValidateMOI', GNOMAD_DOM_HOM_THRESHOLD])
        self.hemi_threshold = config_retrieve(['ValidateMOI', GNOMAD_HEMI_THRESHOLD])

        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS}
            | {key: self.hemi_threshold for key in INFO_HEMI}
            | {'gnomad_ac': self.ac_threshold, 'gnomad_af': self.ad_threshold},
            StructuralVariant.__name__: {key: self.hom_threshold for key in SV_HOMS}
            | {key: self.hemi_threshold for key in SV_HEMI},
        }

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        if variant is present and sufficiently rare, we take it
        discarded if support

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # never apply dominant MOI to support variants
        # more stringent Pop.Freq checks for dominant - hemi restriction
        if (
            principal.support_only
            or self.check_frequency_fails(
                principal.info,
                self.freq_tests[principal.__class__.__name__],
            )
            or self.check_callset_af_fails(principal.info, self.ad_threshold)
        ):
            return classifications

        # all samples which have a variant call
        samples_with_this_variant = principal.het_samples.union(principal.hom_samples)

        for sample_id in samples_with_this_variant:
            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised
            # force minimum depth
            if (
                not (
                    principal.sample_category_check(sample_id, allow_support=False)
                    and self.pedigree.by_id[sample_id].affected == '2'
                )
            ) or principal.check_read_depth(sample_id, self.minimum_depth, principal.info.get(CB1)):
                continue

            # check if this is a candidate for dominant inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_with_this_variant,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree.by_id[sample_id].family,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
                ),
            )
        return classifications


class XPseudoDominantFemale(BaseMoi):
    """
    X-Dominant method which only evaluates females, on the basis that a healthy allele could be inactivated
    A special case of X-Dominant, where we consider Het. variants, used to assess genes which are on X and in genes
    associated with a Biallelic MOI

    Basically a Dominant MOI to be applied to Recessive genes, and results will be labelled as cautionary
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'X_PseudoDominant'):
        """
        accept male hets and homs, and female hets without support

        Args:
            pedigree ():
            applied_moi ():
        """
        self.ad_threshold = config_retrieve(['ValidateMOI', GNOMAD_RARE_THRESHOLD])
        self.ac_threshold = config_retrieve(['ValidateMOI', GNOMAD_AD_AC_THRESHOLD])
        self.hom_threshold = config_retrieve(['ValidateMOI', GNOMAD_DOM_HOM_THRESHOLD])

        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_threshold for key in INFO_HOMS}
            | {'gnomad_ac': self.ac_threshold, 'gnomad_af': self.ad_threshold},
            StructuralVariant.__name__: {key: self.hom_threshold for key in SV_HOMS},
        }

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        if variant is present and sufficiently rare, we take it
        discarded if support

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """
        # unused in this class, we always run this with partial penetrance
        _unused = partial_pen

        classifications = []

        if principal.support_only:
            return classifications

        # never apply dominant MOI to support variants
        # more stringent Pop.Freq checks for dominant - hemi restriction
        if self.check_frequency_fails(
            principal.info,
            self.freq_tests[principal.__class__.__name__],
        ):
            return classifications

        # all females which have a variant call
        females_under_consideration = {sam for sam in principal.het_samples if self.pedigree.by_id[sam].sex == '2'}
        all_with_variant = principal.het_samples.union(principal.hom_samples)
        for sample_id in females_under_consideration:
            # skip primary analysis for unaffected members
            # we require this specific sample to be categorised
            # force minimum depth
            if (
                not (
                    principal.sample_category_check(sample_id, allow_support=False)
                    and self.pedigree.by_id[sample_id].affected == '2'
                )
            ) or principal.check_read_depth(sample_id, self.minimum_depth, principal.info.get(CB1)):
                continue

            # check if this is a candidate for dominant inheritance
            # as we're allowing for flexible 'penetrance' in females, we send all het and hom variants, but allow for a
            # partial penetrance check - the participants can have the variant, but not be affected. They may not be
            # affected without the variant call.
            # There's a slight breakdown here as the males should be interpreted under a full penetrance model, and
            # females under partial penetrance, but that's not trivial without creating a second familial check method.
            # Leaving that aside now as the current implementation is pretty central to the algorithm. Will revisit if
            # this is noisy.
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=all_with_variant,
                partial_pen=True,
            ):
                continue

            classifications.append(
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree.by_id[sample_id].family,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    reasons={self.applied_moi},
                    genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                    flags=principal.get_sample_flags(sample_id)
                    | {'Affected female with heterozygous variant in XLR gene'},
                    independent=True,
                ),
            )
        return classifications


class XRecessiveMale(BaseMoi):
    """
    accept all male non-ref GTs - male variants HOM because GATK
    effectively the same as AutosomalDominant?
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'X_Male'):
        """
        set parameters specific to male X tests

        Args:
            pedigree ():
            applied_moi ():
        """

        self.hom_dom_threshold = config_retrieve(['ValidateMOI', GNOMAD_DOM_HOM_THRESHOLD])
        self.hemi_threshold = config_retrieve(['ValidateMOI', GNOMAD_HEMI_THRESHOLD])

        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_dom_threshold for key in INFO_HOMS}
            | {key: self.hemi_threshold for key in INFO_HEMI},
            StructuralVariant.__name__: {key: self.hom_dom_threshold for key in SV_HOMS}
            | {key: self.hemi_threshold for key in SV_HEMI},
        }

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """
        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if self.check_frequency_fails(
            principal.info,
            self.freq_tests[principal.__class__.__name__],
        ):
            return classifications

        # combine het and hom here, we don't trust the variant callers
        # if hemi count is too high, don't consider males
        # never consider support variants on X for males
        males = {
            sam for sam in principal.het_samples.union(principal.hom_samples) if self.pedigree.by_id[sam].sex == '1'
        }

        for sample_id in males:
            # specific affected sample category check, never consider support on X for males
            if (
                not (
                    self.pedigree.by_id[sample_id].affected == '2'
                    and principal.sample_category_check(sample_id, allow_support=False)
                )
            ) or principal.check_read_depth(sample_id, self.minimum_depth, principal.info.get(CB1)):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(sample_id=sample_id, called_variants=males, partial_pen=partial_pen):
                continue

            classifications.append(
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree.by_id[sample_id].family,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
                ),
            )
        return classifications


class XRecessiveFemaleHom(BaseMoi):
    """
    only consider HOM females
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'X_Recessive HOM Female'):
        """
        set parameters specific to recessive tests

        Args:
            pedigree ():
            applied_moi ():
        """

        self.hom_rec_threshold = config_retrieve(['ValidateMOI', GNOMAD_REC_HOM_THRESHOLD])
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_rec_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {key: self.hom_rec_threshold for key in SV_HOMS},
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,  # noqa: ARG002
        partial_pen: bool = False,
    ) -> list[ReportVariant]:
        """

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        classifications = []

        # remove from analysis if too many homs are present in population databases
        if principal.support_only or self.check_frequency_fails(
            principal.info,
            self.freq_tests[principal.__class__.__name__],
        ):
            return classifications

        # never consider support homs
        samples_to_check = {sam for sam in principal.hom_samples if self.pedigree.by_id[sam].sex == '2'}

        for sample_id in samples_to_check:
            # specific affected sample category check
            if (
                not (
                    self.pedigree.by_id[sample_id].affected == '2'
                    and principal.sample_category_check(sample_id, allow_support=False)
                )
            ) or principal.check_read_depth(sample_id, self.minimum_depth, principal.info.get(CB1)):
                continue

            # check if this is a possible candidate for homozygous inheritance
            if not self.check_familial_inheritance(
                sample_id=sample_id,
                called_variants=samples_to_check,
                partial_pen=partial_pen,
            ):
                continue

            classifications.append(
                ReportVariant(
                    sample=sample_id,
                    family=self.pedigree.by_id[sample_id].family,
                    gene=principal.info.get('gene_id'),
                    var_data=principal,
                    categories=principal.category_values(sample_id),
                    genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                    reasons={self.applied_moi},
                    flags=principal.get_sample_flags(sample_id),
                    independent=True,
                ),
            )
        return classifications


class XRecessiveFemaleCH(BaseMoi):
    """
    ignore males, accept female comp-het only
    """

    def __init__(self, pedigree: Pedigree, applied_moi: str = 'X_RecessiveFemaleCompHet'):
        """
        set parameters specific to recessive tests

        Args:
            pedigree ():
            applied_moi ():
        """

        self.hom_rec_threshold = config_retrieve(['ValidateMOI', GNOMAD_REC_HOM_THRESHOLD])
        self.freq_tests = {
            SmallVariant.__name__: {key: self.hom_rec_threshold for key in INFO_HOMS},
            StructuralVariant.__name__: {key: self.hom_rec_threshold for key in SV_HOMS},
        }
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: VARIANT_MODELS,
        comp_het: CompHetDict | None = None,
        partial_pen: bool = False,  # noqa: ARG002
    ) -> list[ReportVariant]:
        """

        Args:
            principal ():
            comp_het ():
            partial_pen ():
        """

        if comp_het is None:
            comp_het = {}
        classifications = []

        # remove from analysis if too many homs are present in population databases
        if self.check_frequency_fails(
            principal.info,
            self.freq_tests[principal.__class__.__name__],
        ):
            return classifications
        het_females = {sam for sam in principal.het_samples if self.pedigree.by_id[sam].sex == '2'}

        # if het females are present, try and find support
        for sample_id in het_females:
            # don't run primary analysis for unaffected
            # we require this specific sample to be categorised - check Cat 4 contents
            if (
                not (
                    self.pedigree.by_id[sample_id].affected == '2'
                    and principal.sample_category_check(sample_id, allow_support=True)
                )
            ) or principal.check_read_depth(sample_id, self.minimum_depth, principal.info.get(CB1)):
                continue

            for partner in check_for_second_hit(
                first_variant=principal.coordinates.string_format,
                comp_hets=comp_het,
                sample=sample_id,
                require_non_support=principal.sample_support_only(sample_id),
            ):
                # allow for de novo check - also screen out high-AF partners
                if not partner.sample_category_check(sample_id, allow_support=True) or self.check_frequency_fails(
                    partner.info,
                    self.freq_tests[partner.__class__.__name__],
                ):
                    continue

                # check for minimum depth in partner
                if partner.check_read_depth(sample_id, self.minimum_depth, partner.info.get(CB1)):
                    continue

                if not self.check_comp_het(sample_id=sample_id, variant_1=principal, variant_2=partner):
                    continue

                classifications.append(
                    ReportVariant(
                        sample=sample_id,
                        family=self.pedigree.by_id[sample_id].family,
                        gene=principal.info.get('gene_id'),
                        var_data=principal,
                        categories=principal.category_values(sample_id),
                        reasons={self.applied_moi},
                        genotypes=self.get_family_genotypes(variant=principal, sample_id=sample_id),
                        # needs to comply with Seqr
                        support_vars={partner.info['seqr_link']},
                        flags=principal.get_sample_flags(sample_id) | partner.get_sample_flags(sample_id),
                        independent=False,
                    ),
                )

        return classifications
