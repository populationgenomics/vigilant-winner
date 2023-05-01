"""
creates a data model for generation of realistic variant data

The data produced by this method can be slotted into an AIP run
(or any other process using a VEP annotated MT) to generate realistic
variant data.

A simple use case would be replacing the current 'semisynthetic' data
in AIP unit tests with a more full-featured batch of test data.

'Semisynthetic' data in use so far has been generated by taking real
annotated data run through AIP, taking the VCF export, and then doing
some obfuscation (changing variants, genotypes, annotations) without
altering the valid data format. This is realistic, but manual and will
not scale well.

A slightly different use case would be an end-to-end test of the
AIP pipeline, where the input data is generated by this script,
and the output data is compared to the expected output data.

Expected results in this case would be determined by the gene ID(s),
the annotation(s), the genotype(s), and the sample affection status.
"""

# pylint: disable=invalid-name,too-many-instance-attributes


import json
from os.path import join
from dataclasses import dataclass, field

import hail as hl
from cloudpathlib import AnyPath

from reanalysis.utils import CustomEncoder

schema = (
    'struct{'
    'locus:str,alleles:array<str>,rsid:str,qual:float64,filters:set<str>,'
    'AC:int32,AF:float64,AN:int32,info:struct{AC:int32,AF:float64,AN:int32},'
    'gnomad_genomes:struct{AF:float64,AN:int32,AC:int32,Hom:int32,Hemi:int32},'
    'gnomad_exomes:struct{AF:float64,AN:int32,AC:int32,Hom:int32,Hemi:int32},'
    'splice_ai:struct{delta_score:float32,splice_consequence:str},'
    'cadd:struct{PHRED:float32},'
    'dbnsfp:struct{REVEL_score:str,MutationTaster_pred:str},'
    'clinvar:struct{clinical_significance:str,gold_stars:int32,allele_id:int32},'
    'vep:struct{'
    'transcript_consequences:array<struct{gene_symbol:str,gene_id:str,'
    'variant_allele:str,consequence_terms:array<str>,transcript_id:str,'
    'protein_id:str,gene_symbol_source:str,canonical:int32,cdna_start:int32,'
    'cds_start:int32,cds_end:int32,biotype:str,protein_start:int32,protein_end:int32,'
    'sift_score:float64,sift_prediction:str,polyphen_prediction:str,polyphen_score:'
    'float64,mane_select:str,lof:str}>,variant_class:str},geneIds:set<str>'
    '}'
)


@dataclass
class BaseFields:
    """
    base row fields
    """

    def __init__(
        self,
        locus: str,
        alleles: list[str],
        filters: set[str] | None = None,
        rsid: str | None = '.',
        qual: float | None = 60.0,
        ac: int | None = 1,
        af: float | None = 0.001,
        an: int | None = 1,
    ):
        self.locus = locus
        self.alleles = alleles
        self.filters = filters or set()
        self.rsid = rsid
        self.qual = qual
        # shouldn't be required in info anymore
        self.info = {'AC': ac, 'AF': af, 'AN': an}
        self.AC = ac
        self.AF = af
        self.AN = an


@dataclass
class AFGeneric:
    """
    generic Allele Frequency data model
    """

    AF: float = field(default=0.001)
    AN: int = field(default=1)
    AC: int = field(default=1)
    Hom: int = field(default=0)


@dataclass
class AFData:
    """
    specific Allele Frequency data model
    """

    gnomad_genomes: AFGeneric = field(default_factory=AFGeneric)
    gnomad_exomes: AFGeneric = field(default_factory=AFGeneric)


@dataclass
class Splice:
    """
    Splice data model
    """

    delta_score: float = field(default=0.01)
    splice_consequence: str = field(default='none')


@dataclass
class CADD:
    """
    CADD data model
    """

    PHRED: float = field(default=0.01)


@dataclass
class DBnsfp:
    """
    DBnsfp data model
    """

    REVEL_score: str = field(default='0.0')
    MutationTaster_pred: str = field(default='n')


@dataclass
class Clinvar:
    """
    Clinvar data model
    """

    clinical_significance: str = field(default_factory=str)
    gold_stars: int = field(default_factory=int)
    allele_id: int = field(default_factory=int)


@dataclass
class TXFields:
    """
    TX fields data model
    """

    gene_symbol: str
    gene_id: str
    variant_allele: str = field(default_factory=str)
    consequence_terms: list = field(default_factory=list)
    transcript_id: str = field(default_factory=str)
    protein_id: str = field(default_factory=str)
    gene_symbol_source: str = field(default_factory=str)
    canonical: int = field(default=1)
    cdna_start: int = field(default=1)
    cds_end: int = field(default=1)
    biotype: str = field(default_factory=str)
    protein_start: int = field(default=1)
    protein_end: int = field(default=1)
    sift_score: int = field(default=1.0)  # lowest possible score
    sift_prediction: str = field(default_factory=str)
    polyphen_prediction: str = field(default='neutral')
    polyphen_score: float = field(default=0.01)
    mane_select: str = field(default_factory=str)
    lof: str = field(default_factory=str)


@dataclass
class Entry:
    """
    entry data model
    """

    def __init__(
        self,
        gt: str,
        ad: list[int] | None = None,
        gq: int | None = None,
        pl: list[int] | None = None,
        ps: int | None = None,
    ):
        """
        This is the per-call data

        Args:
            gt (str): needs to be parsed into a call later, not JSON'able
            ad (list[int]): depths for normalised alleles
            gq (int): genotype quality
            pl (list[int]): phred-scaled likelihoods
                auto-selected based on GQ & GT if not provided
        """
        self.GT = gt
        self.AD = ad or [15, 15]
        self.DP = sum(self.AD)
        self.GQ = gq or 60
        self.PL = (
            pl
            or [[0, self.GQ, 1000], [self.GQ, 0, 1000], [1000, self.GQ, 0]][
                gt.count('1')
            ]
        )
        self.PS = ps

    @staticmethod
    def get_schema_entry():
        """
        how to represent this data type
        Returns:

        """
        return (
            'struct{GT:str,AD:array<int32>,DP:int32,GQ:int32,PL:array<int32>,PS:int32}'
        )


@dataclass
class VepVariant:
    """
    class object to sweep up all the data models
    """

    def __init__(
        self,
        base: BaseFields,
        tx: list[TXFields],
        sample_data: dict[str, Entry] | None = None,
        af: AFData | None = None,
        cadd: CADD | None = None,
        dbnsfp: DBnsfp | None = None,
        clinvar: Clinvar | None = None,
        splice: Splice | None = None,
        var_class: str = 'SNV',
    ):
        """
        VEP variant data model
        Args:
            tx (list[TXFields]): one or more transcript consequences
            base (BaseFields): global variables
            sample_data (dict[str, Entry]): per-sample data
            af (AFData): Allele Freq data, or None
            cadd (CADD): CADD data, or None
            dbnsfp (DBnsfp): Revel/MutationTaster, or None
            clinvar (Clinvar): Clinvar data, or None
            splice (Splice): SpliceAI data, or None
            var_class (str): ... SNV
        """

        self.data = (
            {
                'vep': {'transcript_consequences': tx, 'variant_class': var_class},
                'geneIds': {tx.gene_id for tx in tx},
                'dbnsfp': dbnsfp or DBnsfp(),
                'cadd': cadd or CADD(),
                'clinvar': clinvar or Clinvar(),
                'splice_ai': splice or Splice(),
            }
            | base.__dict__
            | (af or AFData()).__dict__
        )

        if sample_data:
            self.data |= sample_data

    def to_string(self) -> str:
        """
        convert this object to a string
        Returns:
            a single line JSON string
        """
        return json.dumps(self.data, cls=CustomEncoder) + '\n'


class SneakyTable:
    """
    class to take multiple individual variants
    and generate a Hail Matrix Table from them
    """

    def __init__(
        self,
        variants: list[VepVariant],
        tmp_path: str,
        sample_details: dict[str, str] | None = None,
    ):
        """
        Args:
            variants (list[VepVariant]): list of VepVariant objects
            sample_details (dict[str,str]): sample IDs and corresponding schema
            tmp_path (str): where to write the hail table
        """
        self.variants = variants
        self.sample_details = sample_details
        self.tmp_path = tmp_path
        try:
            hl.init(default_reference='GRCh38')
        except BaseException:  # pylint: disable=broad-except
            pass

    def modify_schema(self) -> str:
        """
        inserts the per-sample data into the schema

        Returns:
            updated schema
        """

        if not self.sample_details:
            return schema

        new_schema = schema[:-1]
        for sample, dtype in self.sample_details.items():
            new_schema += f',{sample}:{dtype}'
        new_schema += '}'
        return new_schema

    def json_to_file(self) -> str:
        """
        write the data to a temp file
        Returns:
            str: path to the temp file
        """

        # write all variants in this object
        json_temp = join(self.tmp_path, 'vep.json')
        with AnyPath(json_temp).open('w') as f:
            for variant in self.variants:
                f.write(variant.to_string())
        return json_temp

    def to_hail(self, hail_table: bool = False) -> hl.MatrixTable | hl.Table:
        """
        write the data model to a hail table

        entertain 3 behaviours:
         - write to a MatrixTable, including the genotypes/entries
         - write to a MatrixTable, but only include the annotations
         - keep as a Table (containing samples or not)

         Args:
            hail_table (bool): if True, return a Table, else a MatrixTable

        Returns:
            the MatrixTable of the faux annotation data
        """

        # update the schema if sample entries were added
        sample_schema = self.modify_schema()
        json_schema = hl.dtype(sample_schema)

        # read JSON data from a hail table
        # field must be f0 if no header
        ht = hl.import_table(
            self.json_to_file(), no_header=True, types={'f0': json_schema}
        )

        # unwrap the data
        ht = ht.transmute(**ht.f0)

        # transmute the locus and alleles, set as keys
        ht = ht.transmute(locus=hl.parse_locus(ht.locus), alleles=ht.alleles).key_by(
            'locus', 'alleles'
        )

        # stopping point for table-only
        if hail_table:
            # checkpoint out to a temp path
            tmp_ht = join(self.tmp_path, 'vep.ht')
            ht.write(tmp_ht, overwrite=True)
            return ht

        if not self.sample_details:
            return hl.MatrixTable.from_rows_table(ht)

        tmp_mt = join(self.tmp_path, 'vep.mt')

        # convert to a matrix table, with sample IDs as columns
        mt = ht.to_matrix_table_row_major(
            columns=list(self.sample_details.keys()), col_field_name='s'
        )

        # parse the genotype calls as hl.call
        mt = mt.annotate_entries(GT=hl.parse_call(mt.GT))
        mt.write(tmp_mt, overwrite=True)

        # send it
        return mt
