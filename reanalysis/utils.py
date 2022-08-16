"""
a collection of classes and methods
which may be shared across reanalysis components
"""

from collections import defaultdict
from dataclasses import dataclass, is_dataclass
from enum import Enum
from itertools import combinations_with_replacement
from typing import Any, Union

import json
import logging
import re

from cloudpathlib import AnyPath


InfoDict = dict[str, Union[str, dict[str, str]]]
PanelAppDict = dict[str, dict[str, Union[str, bool]]]

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: set[int] = {HOMREF, UNKNOWN}

PHASE_SET_DEFAULT = -2147483648
X_CHROMOSOME = {'X'}
NON_HOM_CHROM = {'Y', 'MT', 'M'}


class FileTypes(Enum):
    """
    enumeration of permitted input file types
    """

    HAIL_TABLE = '.ht'
    MATRIX_TABLE = '.mt'
    VCF = '.vcf'
    VCF_GZ = '.vcf.gz'
    VCF_BGZ = '.vcf.bgz'


def identify_file_type(file_path: str) -> FileTypes | Exception:
    """
    return type of the file, if present in FileTypes enum

    :param file_path:
    :return:
    """
    pl_filepath = AnyPath(file_path)

    assert (
        pl_filepath.exists()
    ), f'input file {file_path} does not exist or is not visible from this host'

    # pull all extensions (e.g. .vcf.bgz will be split into [.vcf, .bgz]
    extensions = pl_filepath.suffixes

    assert len(extensions) > 0, 'cannot identify input type from extensions'

    if extensions[-1] == '.ht':
        return FileTypes.HAIL_TABLE
    if extensions[-1] == '.mt':
        return FileTypes.MATRIX_TABLE
    if extensions == ['.vcf']:
        return FileTypes.VCF
    if extensions == ['.vcf', '.gz']:
        return FileTypes.VCF_GZ
    if extensions == ['.vcf', '.bgz']:
        return FileTypes.VCF_BGZ
    raise Exception(f'File cannot be definitively typed: {str(extensions)}')


@dataclass
class Coordinates:
    """
    a home for the positional variant attributes
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def string_format(self):
        """
        forms a string representation
        chr-pos-ref-alt
        """
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'


@dataclass
class AbstractVariant:  # pylint: disable=too-many-instance-attributes
    """
    create a bespoke variant class
    pull all content out of the cyvcf2 object
    we could have a separate implementation for pyvcf,
    or for a direct parser...
    """

    def __init__(
        self,
        var,
        samples: list[str],
        config: dict[str, Any],
        as_singletons=False,
    ):
        """
        var is a cyvcf2.Variant
        :param var:
        :param samples:
        :param config:
        :param as_singletons:
        """

        # extract the coordinates into a separate object
        self.coords = Coordinates(
            var.CHROM.replace('chr', ''), var.POS, var.REF, var.ALT[0]
        )

        # set the class attributes
        self.category_1: bool = var.INFO.get('Category1') == 1
        self.category_2: bool = var.INFO.get('Category2') == 1
        self.category_3: bool = var.INFO.get('Category3') == 1
        self.category_5: bool = var.INFO.get('Category5') == 1

        # de novo categories are a list of strings or empty list
        # if cohort runs as singletons, remove possibility of de novo
        if as_singletons:
            self.category_4 = []
            self.category_4b = []
        else:
            self.category_4: list[str] = (
                var.INFO.get('Category4').split(',')
                if var.INFO.get('Category4') != 'missing'
                else []
            )
            self.category_4b: list[str] = (
                var.INFO.get('Category4b').split(',')
                if var.INFO.get('Category4b') != 'missing'
                else []
            )

        self.category_support: bool = var.INFO.get('CategorySupport') == 1

        # get all zygosities once per variant
        # abstraction avoids pulling per-sample calls again later
        self.het_samples, self.hom_samples = get_non_ref_samples(
            variant=var, samples=samples
        )

        self.info: dict[str, str] = extract_info(variant=var)
        self.transcript_consequences: list[dict[str, str]] = extract_csq(
            variant=var, config=config
        )

        # identify variant sets phased with this one
        # cyvcf2 uses a default value for the phase set, skip that
        # this is restricted to a single int for phase_set
        self.phased = defaultdict(dict)

        # first set the numpy.ndarray to be a list of ints
        # the zip against ordered sample IDs
        # this might need to store the exact genotype too
        # i.e. 0|1 and 1|0 can be in the same phase-set
        # but are un-phased variants
        for sample, phase, genotype in zip(
            samples, map(int, var.format('PS')), var.genotypes
        ):
            # cyvcf2.Variant holds two ints, and a bool
            allele_1, allele_2, phased = genotype
            if not phased:
                continue
            gt = f'{allele_1}|{allele_2}'
            # phase set is a number
            if phase != PHASE_SET_DEFAULT:
                self.phased[sample][phase] = gt

        self.ab_ratios = dict(zip(samples, map(float, var.gt_alt_freqs)))

    def __str__(self):
        return repr(self)

    @property
    def category_1_2_3_5(self) -> bool:
        """
        check that the variant has at least one assigned class
        :return:
        """
        return any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
                self.category_5,
            ]
        )

    @property
    def category_non_support(self) -> bool:
        """
        check the variant has at least one non-support category assigned
        :return:
        """
        return any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
                self.category_4,
                self.category_4b,
                self.category_5,
            ]
        )

    @property
    def is_classified(self) -> bool:
        """
        check that the variant has at least one assigned class
        supporting category is considered here
        :return:
        """
        return any(
            [
                self.category_1,
                self.category_2,
                self.category_3,
                self.category_4,
                self.category_4b,
                self.category_5,
                self.category_support,
            ]
        )

    @property
    def support_only(self) -> bool:
        """
        checks that the variant was only class 4
        :return:
        """
        return self.category_support and not self.category_non_support

    def category_values(self, sample: str) -> list[str]:
        """
        get a list values representing the classes present on this variant
        for each category, append that number if the class is present
        - support is not an int
        - de novo on a per-sample basis
        """

        # for basic categories, append the value if flag is present
        categories = [
            value
            for value, category in [
                ('1', self.category_1),
                ('2', self.category_2),
                ('3', self.category_3),
                ('5', self.category_5),
                ('in_silico', self.category_support),
            ]
            if category
        ]

        # sample-specific check for de novo
        if self.sample_de_novo(sample_id=sample):
            categories.append('de_novo')

        return categories

    def sample_de_novo(self, sample_id: str) -> bool:
        """
        takes a specific sample ID, to check if the sample has a de novo call

        :param sample_id:
        :return:
        """

        return sample_id in self.category_4 or sample_id in self.category_4b

    def sample_specific_category_check(self, sample_id: str) -> bool:
        """

        :param sample_id:
        :return:
        """
        return self.category_1_2_3_5 or self.sample_de_novo(sample_id)

    def get_sample_flags(self, sample: str) -> list[str]:
        """
        gets all report flags for this sample
        """
        return self.check_ab_ratio(sample) + self.check_dodgy_de_novo(sample)

    def check_dodgy_de_novo(self, sample: str) -> list[str]:
        """
        flag if a de novo is only called by the lenient method
        """
        if self.sample_de_novo(sample) and sample not in self.category_4b:
            return ['Dodgy de novo']
        return []

    def check_ab_ratio(self, sample: str) -> list[str]:
        """
        AB ratio test for this sample's variant call

        ab = matrix.AD[1] / hl.sum(matrix.AD)
        return matrix.filter_entries(
            (matrix.GT.is_hom_ref() & (ab <= 0.15))
            | (matrix.GT.is_het() & (ab >= 0.25) & (ab <= 0.75))
            | (matrix.GT.is_hom_var() & (ab >= 0.85))
        )
        :param sample: this affected individual
        """
        het = sample in self.het_samples
        hom = sample in self.hom_samples
        variant_ab = self.ab_ratios.get(sample, 0.0)
        if (
            (variant_ab <= 0.15)
            or (het and not 0.25 <= variant_ab <= 0.75)
            or (hom and variant_ab <= 0.85)
        ):
            return ['AB Ratio']
        return []


# CompHetDict structure: {sample: {variant_string: [variant, ...]}}
# sample: string, e,g, CGP12345
CompHetDict = dict[str, dict[str, list[AbstractVariant]]]
GeneDict = dict[str, list[AbstractVariant]]


@dataclass
class ReportedVariant:
    """
    minimal model representing variant categorisation event
    the initial variant
    the MOI passed
    the support (if any)
    allows for the presence of flags e.g. Borderline AB ratio
    """

    sample: str
    gene: str
    var_data: AbstractVariant
    reasons: set[str]
    supported: bool
    support_vars: list[str] | None = None
    flags: list[str] | None = None


def canonical_contigs_from_vcf(reader) -> set[str]:
    """
    reader is a cyvcf2.VCFReader
    read the header fields from the VCF handle
    return a set of all 'canonical' contigs
    :param reader:
    :return:
    """

    # contig matching regex - remove all HLA/decoy/unknown
    contig_re = re.compile(r'^(chr)?[0-9XYMT]{1,2}$')

    return {
        contig['ID']
        for contig in reader.header_iter()
        if contig['HeaderType'] == 'CONTIG' and re.match(contig_re, contig['ID'])
    }


def gather_gene_dict_from_contig(
    contig: str,
    variant_source,
    config: dict[str, Any],
    panelapp_data: PanelAppDict,
    singletons: bool = False,
    blacklist: list[str] | None = None,
) -> GeneDict:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup
    {
        gene1: [var1, var2],
        gene2: [var3],
        ...
    }
    :param contig: contig name from header (canonical_contigs_from_vcf)
    :param variant_source: VCF reader instance
    :param config: configuration file
    :param panelapp_data:
    :param singletons:
    :param blacklist:
    :return: populated lookup dict
    """

    if blacklist is None:
        blacklist = []

    # a dict to allow lookup of variants on this whole chromosome
    contig_variants = 0
    contig_dict = defaultdict(list)

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant in variant_source(contig):
        abs_var = AbstractVariant(
            var=variant,
            samples=variant_source.samples,
            config=config,
            as_singletons=singletons,
        )

        if abs_var.coords.string_format in blacklist:
            logging.info(
                f'Skipping blacklisted variant: {abs_var.coords.string_format}'
            )
            continue

        # if the gene isn't 'new' in PanelApp, remove Class2 flag
        if abs_var.category_2:
            gene_id = abs_var.info.get('gene_id')
            gene_data = panelapp_data.get(gene_id, False)
            if not gene_data:
                continue

            if not gene_data.get('new', False):
                variant.category_2 = False

        # if unclassified, skip the whole variant
        if not abs_var.is_classified:
            continue

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict[abs_var.info.get('gene_id')].append(abs_var)

    logging.info(f'Contig {contig} contained {contig_variants} variants')
    logging.info(f'Contig {contig} contained {len(contig_dict)} genes')
    return contig_dict


def read_json_from_path(bucket_path: str) -> dict[str, Any]:
    """
    take a path to a JSON file, read into an object
    :param bucket_path:
    """
    with AnyPath(bucket_path).open() as handle:
        return json.load(handle)


def get_simple_moi(panel_app_moi: str) -> str:
    """
    takes the vast range of PanelApp MOIs, and reduces to a reduced
    range of cases which can be easily implemented in RD analysis
    This is required to reduce the complexity of an MVP
    Could become a strict enumeration

    {
        Biallelic: [all, panelapp, moi, equal, to, biallelic],
        Monoallelic: [all MOI to be interpreted as monoallelic]
    }

    :param panel_app_moi: full PanelApp string
    :return: a simplified representation
    """

    # default to considering both. NOTE! Many genes have Unknown MOI!
    simple_moi = 'Mono_And_Biallelic'
    if panel_app_moi is None or panel_app_moi == 'Unknown':
        # exit iteration, all simple moi considered
        return simple_moi

    # ideal for match-case, coming to a python 3.10 near you!
    if panel_app_moi.startswith('BIALLELIC'):
        simple_moi = 'Biallelic'
    if panel_app_moi.startswith('BOTH'):
        simple_moi = 'Mono_And_Biallelic'
    if panel_app_moi.startswith('MONO'):
        simple_moi = 'Monoallelic'
    if panel_app_moi.startswith('X-LINKED'):
        if 'biallelic' in panel_app_moi:
            simple_moi = 'Hemi_Bi_In_Female'
        else:
            simple_moi = 'Hemi_Mono_In_Female'

    return simple_moi


def get_non_ref_samples(variant, samples: list[str]) -> tuple[set[str], set[str]]:
    """
    variant is a cyvcf2.Variant
    for this variant, find all samples with a call
    cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    return het, hom, and the union of het and hom
    :param variant:
    :param samples:
    :return:
    """
    het_samples = set()
    hom_samples = set()

    # this iteration is based on the cyvcf2 representations
    for sam, genotype_int in zip(samples, variant.gt_types):

        if genotype_int in BAD_GENOTYPES:
            continue
        if genotype_int == 1:
            het_samples.add(sam)
        if genotype_int == 3:
            hom_samples.add(sam)

    return het_samples, hom_samples


def extract_csq(variant, config: dict[str, dict[str, str]]):
    """
    variant is a cyvcf2.Variant
    specifically handle extraction of the CSQ list
    :param variant:
    :param config:
    :return:
    """

    # pull the CSQ content from the variant
    csq_contents = variant.INFO.get('CSQ')

    # allow for no CSQ data, i.e. splice variant
    if not csq_contents:
        return []

    # config region concerning variant objects
    # break mono-CSQ-string into components
    csq_categories = config['variant_object']['csq_string']

    # iterate over all consequences, and make each into a dict
    return [
        dict(zip(csq_categories, each_csq.split('|')))
        for each_csq in csq_contents.split(',')
    ]


def extract_info(variant):
    """
    variant is a cyvcf2.Variant
    creates an INFO dict by pulling content from the variant info
    keeps a list of dictionaries for each transcript_consequence
    :param variant:
    :return:
    """

    # choose some values to exclude, and keep everything else
    exclusions = {
        'csq',
        'category1',
        'category2',
        'category3',
        'category4',
        'category4b',
        'categorysupport',
        'support_only',
    }

    # grab the basic information from INFO
    return {x.lower(): y for x, y in variant.INFO if x.lower() not in exclusions}


class CustomEncoder(json.JSONEncoder):
    """
    to be used as a JSON encoding class
    - replaces all sets with lists
    - replaces dataclass objects with a dictionary of the same
    """

    def default(self, o):
        """
        takes an arbitrary object, and forms a JSON representation
        where the object doesn't have an easy string representation,
        transform to a valid object: set -> list, class -> dict
        :param o:
        :return:
        """

        if is_dataclass(o):
            return o.__dict__
        if isinstance(o, set):
            return list(o)
        return json.JSONEncoder.default(self, o)


def find_comp_hets(var_list: list[AbstractVariant], pedigree) -> CompHetDict:
    """
    manual implementation to find compound hets
    variants provided in the format

    [var1, var2, ..]

    generate pair content in the form
    {
        sample: {
            var_as_string: [partner_variant, ...]
        }
    }

    :param var_list:
    :param pedigree: peddy.Ped
    :return:
    """

    # create an empty dictionary
    comp_het_results = defaultdict(dict)

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coords.chrom == var_2.coords.chrom

        if (var_1.coords == var_2.coords) or var_1.coords.chrom in NON_HOM_CHROM:
            continue

        sex_chrom = var_1.coords.chrom in X_CHROMOSOME

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False
            ped_sample = pedigree.get(sample)

            # don't assess male compound hets on sex chromosomes
            if ped_sample.sex == 'male' and sex_chrom:
                continue

            # check for both variants being in the same phase set
            if sample in var_1.phased and sample in var_2.phased:

                # check for presence of the same phase set
                for phase_set in [
                    ps
                    for ps in var_1.phased[sample].keys()
                    if ps in var_2.phased[sample].keys()
                ]:
                    if (
                        var_1.phased[sample][phase_set]
                        == var_2.phased[sample][phase_set]
                    ):
                        phased = True
            if not phased:
                comp_het_results[sample].setdefault(
                    var_1.coords.string_format, []
                ).append(var_2)
                comp_het_results[sample].setdefault(
                    var_2.coords.string_format, []
                ).append(var_1)

    return comp_het_results
