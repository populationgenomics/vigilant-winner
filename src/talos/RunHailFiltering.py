"""
Read, filter, annotate, classify, and write Genetic data
- read MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC)
- extract generic fields
- remove all rows and consequences not relevant to GREEN genes
- extract vep data into CSQ string(s)
- annotate with categories 1, 2, 3, 4, 5, 6, and Support
- remove un-categorised variants
- write as VCF
"""

from argparse import ArgumentParser

import hail as hl
from peds import open_ped

from talos.config import config_retrieve
from talos.hail_audit import (
    BASE_FIELDS_REQUIRED,
    FIELDS_REQUIRED,
    USELESS_FIELDS,
    VEP_TX_FIELDS_REQUIRED,
    fields_audit,
    vep_audit,
)
from talos.models import PanelApp
from talos.static_values import get_logger
from talos.utils import read_json_from_path

# set some Hail constants
MISSING_INT = hl.int32(0)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_STRING = hl.str('missing')
ONE_INT = hl.int32(1)
BENIGN = hl.str('benign')
LOFTEE_HC = hl.str('HC')
PATHOGENIC = hl.str('pathogenic')
SPLICE_ALTERING = hl.str('splice-altering')

# decide whether to repartition the data before processing starts
MAX_PARTITIONS = 10000


def annotate_clinvarbitration(mt: hl.MatrixTable, clinvar: str) -> hl.MatrixTable:
    """
    Don't allow these annotations to be missing
    - Talos has been co-developed with ClinvArbitration, a ClinVar re-summary effort
    - We no longer permit this to be missing (this has slipped in the past, causing odd results)

    See: https://github.com/populationgenomics/ClinvArbitration

    We replace any existing ClinVar annotations with our own version, then identify Pathogenic/Benign variants

    Args:
        mt (): the MatrixTable of all variants
        clinvar (): the table of private ClinVar annotations to use

    Returns:
        The same MatrixTable but with additional annotations
    """

    get_logger().info(f'loading private clinvar annotations from {clinvar}')
    ht = hl.read_table(clinvar)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_significance=hl.or_else(ht[mt.row_key].clinical_significance, MISSING_STRING),
            clinvar_stars=hl.or_else(ht[mt.row_key].gold_stars, MISSING_INT),
            clinvar_allele=hl.or_else(ht[mt.row_key].allele_id, MISSING_INT),
        ),
    )

    # remove all confidently benign
    mt = mt.filter_rows(
        (mt.info.clinvar_significance.lower().contains(BENIGN)) & (mt.info.clinvar_stars > 0),
        keep=False,
    )

    # annotate as either strong or regular, return the result
    return mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_talos=hl.if_else(
                mt.info.clinvar_significance.lower().contains(PATHOGENIC),
                ONE_INT,
                MISSING_INT,
            ),
            clinvar_talos_strong=hl.if_else(
                (mt.info.clinvar_significance.lower().contains(PATHOGENIC)) & (mt.info.clinvar_stars > 0),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def annotate_exomiser(mt: hl.MatrixTable, exomiser: str | None = None) -> hl.MatrixTable:
    """
    Annotate this MT with top hits from Exomiser

    The exomiser table must be indexed on [locus, alleles], with an extra column "proband_details"

    Args:
        mt (): the MatrixTable of all variants
        exomiser (str | None): optional, path to HT of exomiser results

    Returns:
        The same MatrixTable but with additional annotations
    """

    if not exomiser:
        get_logger().info('No exomiser table found, skipping annotation')
        return mt.annotate_rows(info=mt.info.annotate(categorydetailsexomiser=MISSING_STRING))

    get_logger().info(f'loading exomiser variants from {exomiser}')
    exomiser_ht = hl.read_table(exomiser)
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorydetailsexomiser=hl.or_else(exomiser_ht[mt.row_key].proband_details, MISSING_STRING),
        ),
    )


def annotate_codon_clinvar(mt: hl.MatrixTable, pm5_path: str | None):
    """
    takes the protein indexed clinvar results and matches up against
    the variant data

    this might be a grossly inefficient method...

    The process is:
    1. load up the hail table of all clinvar annotations indexed by residue
    2. re-shuffle the current variant MT to be similarly indexed on residue
    3. match the datasets across to get residues, clinvar, and corresponding
        loci within this callset
    4. explode that back out to index all the clinvar alleles directly on the
        loci relevant to this callset
    5. use that final table to annotate all relevant clinvar allele IDs as a
        flag, indicating when a variant in this callset creates a residue change
        seen in clinvar

    Note - this might include more munging than necessary, but frankly, I'm not
    smart enough to simplify it. I'll bookmark it for future consideration.

    "Everyone knows that debugging is twice as hard as writing a program in the
    first place. So if you're as clever as you can be when you write it, how will
    you ever debug it?” ― Brian Kernighan

    This matching is universal, i.e. if a variant is the exact position and change
    creating a known pathogenic missense, this method should always find that
    annotation through this lookup. It would be problematic if it didn't. For the
    purposes of PM5, it would be invalid to use the exact same allele as further
    evidence, so exact matches must be filtered out downstream.

    Args:
        mt (): MT of all variants
        pm5_path (str or None): path to a (localised) PM5 annotations Hail Table

    Returns:
        Same MT with an extra category label containing links to all clinvar
        missense variants affecting the same residue as a missense in this
        call-set - shared residue affected on at least one transcript
    """

    if pm5_path is None:
        get_logger().info('PM5 table not found, skipping annotation')
        return mt.annotate_rows(info=mt.info.annotate(categorydetailspm5=MISSING_STRING))

    # read in the codon table
    get_logger().info(f'Reading clinvar alleles by codon from {pm5_path}')
    codon_clinvar = hl.read_table(pm5_path)

    # boom those variants out by consequence
    codon_variants = mt.explode_rows(mt.vep.transcript_consequences).rows()

    # filter for missense, no indels (don't trust VEP)
    codon_variants = codon_variants.filter(
        (hl.len(codon_variants.alleles[0]) == ONE_INT)
        & (hl.len(codon_variants.alleles[1]) == ONE_INT)
        & (codon_variants.vep.transcript_consequences.consequence_terms.contains('missense_variant')),
    )

    # set the protein residue as an attribute
    codon_variants = codon_variants.annotate(
        residue_affected=hl.str('::').join(
            [
                codon_variants.vep.transcript_consequences.protein_id,
                hl.str(codon_variants.vep.transcript_consequences.protein_start),
            ],
        ),
    )

    # 4. re-key the table on Transcript::Codon
    codon_variants = codon_variants.key_by(codon_variants.residue_affected)

    # 5. extract the position table (protein change linked to all loci)
    codon_variants = codon_variants.select(codon_variants.locus, codon_variants.alleles).collect_by_key(
        name='positions',
    )

    # join the real variant positions with aggregated clinvar
    # 'values' here is the array of all positions
    codon_variants = codon_variants.join(codon_clinvar)

    # explode back out to release the positions
    codon_variants = codon_variants.explode(codon_variants.positions)

    # annotate positions back to normal names (not required?)
    codon_variants = codon_variants.transmute(
        locus=codon_variants.positions.locus,
        alleles=codon_variants.positions.alleles,
    )

    # re-key by locus/allele
    codon_variants = codon_variants.key_by(codon_variants.locus, codon_variants.alleles)

    # aggregate back to position and alleles
    codon_variants = codon_variants.select(codon_variants.clinvar_alleles).collect_by_key(name='clinvar_variations')

    codon_variants = codon_variants.annotate(
        clinvar_variations=hl.str('+').join(
            hl.set(hl.map(lambda x: x.clinvar_alleles, codon_variants.clinvar_variations)),
        ),
    )

    # conditional annotation back into the original MT
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorydetailspm5=hl.or_else(codon_variants[mt.row_key].clinvar_variations, MISSING_STRING),
        ),
    )


def annotate_splicevardb(mt: hl.MatrixTable, svdb_path: str | None):
    """
    Takes the locus,ref,alt indexed table of SpliceVarDB variants, and matches
    them up against the variant data
    Annotates with a boolean flag if the variant is Splice-altering according to SVDB

    Args:
        mt (): MT of all variants
        svdb_path (str or None): path to a (localised) PM5 annotations Hail Table

    Returns:
        Same MT with an extra category label
    """

    if svdb_path is None:
        get_logger().info('SVDB table not found, skipping annotation')
        return mt.annotate_rows(info=mt.info.annotate(categorydetailsSVDB=MISSING_STRING))

    # read in the codon table
    get_logger().info(f'Reading SpliceVarDB data from {svdb_path}')
    svdb_ht = hl.read_table(svdb_path)

    # annotate relevant variants with the SVDB results
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            svdb_classification=hl.or_else(svdb_ht[mt.row_key].classification, MISSING_STRING),
            svdb_location=hl.or_else(svdb_ht[mt.row_key].location, MISSING_STRING),
            svdb_method=hl.or_else(svdb_ht[mt.row_key].method, MISSING_STRING),
            svdb_doi=hl.or_else(svdb_ht[mt.row_key].doi, MISSING_STRING),
        ),
    )

    # annotate category if Splice-altering according to SVDB
    return mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleansvdb=hl.if_else(
                mt.info.svdb_classification.lower().contains(SPLICE_ALTERING),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def filter_on_quality_flags(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    filter MT to rows with 0 quality filters
    note: in Hail, PASS is represented as an empty set

    This is overridden with Clinvar Pathogenic

    Args:
        mt (hl.MatrixTable): all remaining variants
    Returns:
        MT with all filtered variants removed
    """

    return mt.filter_rows(
        hl.is_missing(mt.filters)
        | (mt.filters.length() == 0)
        | (
            (mt.info.clinvar_talos == ONE_INT)
            | (mt.info.categorybooleansvdb == ONE_INT)
            | (mt.info.categorydetailsexomiser != MISSING_STRING)
        ),
    )


def filter_to_well_normalised(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    single alt per row, no missing Alt

    Args:
        mt (hl.MatrixTable):
    Returns:
        filtered MT
    """
    biallelic = 2
    return mt.filter_rows((hl.len(mt.alleles) == biallelic) & (mt.alleles[1] != '*'))


def extract_annotations(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields which aren't per-consequence
    store these in INFO (required to be included in VCF export)
    replace with placeholder (least consequential) if empty
    e.g. most tools score 0, but for Sift 1 is the least important

    Args:
        mt ():
    Returns:
        Same matrix with re-positioned attributes
    """

    get_logger().info('Pulling VEP annotations into INFO field')

    return mt.annotate_rows(
        info=mt.info.annotate(
            gnomad_ex_af=hl.or_else(mt.gnomad_exomes.AF, MISSING_FLOAT_LO),
            gnomad_ex_an=hl.or_else(mt.gnomad_exomes.AN, MISSING_INT),
            gnomad_ex_ac=hl.or_else(mt.gnomad_exomes.AC, MISSING_INT),
            gnomad_ex_hom=hl.or_else(mt.gnomad_exomes.Hom, MISSING_INT),
            gnomad_ex_hemi=hl.or_else(mt.gnomad_exomes.Hemi, MISSING_INT),
            gnomad_af=hl.or_else(mt.gnomad_genomes.AF, MISSING_FLOAT_LO),
            gnomad_an=hl.or_else(mt.gnomad_genomes.AN, MISSING_INT),
            gnomad_ac=hl.or_else(mt.gnomad_genomes.AC, MISSING_INT),
            gnomad_hom=hl.or_else(mt.gnomad_genomes.Hom, MISSING_INT),
            gnomad_hemi=hl.or_else(mt.gnomad_genomes.Hemi, MISSING_INT),
            splice_ai_delta=hl.or_else(mt.splice_ai.delta_score, MISSING_FLOAT_LO),
            splice_ai_csq=hl.or_else(mt.splice_ai.splice_consequence, MISSING_STRING).replace(' ', '_'),
        ),
    )


def permissive_filter_matrix_by_ac(mt: hl.MatrixTable, ac_threshold: float = 0.01) -> hl.MatrixTable:
    """
    Remove variants with AC in joint-call over threshold
    Will never remove variants with 5 or fewer instances
    Retain high frequency variants if already prioritised in ClinVar or Exomiser

    Args:
        mt (hl.MatrixTable):
        ac_threshold (float):
    Returns:
        MT with all common-in-this-JC variants removed (unless overridden by clinvar path)
    """
    min_callset_ac = 5
    return mt.filter_rows(
        ((min_callset_ac >= mt.info.AC[0]) | (ac_threshold > mt.info.AC[0] / mt.info.AN))
        | ((mt.info.clinvar_talos == ONE_INT) | (mt.info.categorydetailsexomiser != MISSING_STRING)),
    )


def filter_matrix_by_ac(mt: hl.MatrixTable, ac_threshold: float = 0.01) -> hl.MatrixTable:
    """
    Remove variants with AC in joint-call over threshold
    Will never remove variants with 5 or fewer instances
    This is now more strict, not permitting clinvar pathogenics to pass through

    Args:
        mt (hl.MatrixTable):
        ac_threshold (float):
    Returns:
        MT with all common-in-this-JC variants removed
    """
    min_callset_ac = 5
    return mt.filter_rows((min_callset_ac >= mt.info.AC[0]) | (ac_threshold > mt.info.AC[0] / mt.info.AN))


def filter_to_population_rare(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad Exomes and Genomes
    """
    # gnomad exomes and genomes below threshold or missing
    # if missing they were previously replaced with 0.0
    # 'semi-rare' as dominant filters will be more strictly filtered later
    rare_af_threshold = config_retrieve(['RunHailFiltering', 'af_semi_rare'])
    return mt.filter_rows(
        (hl.or_else(mt.gnomad_exomes.AF, MISSING_FLOAT_LO) < rare_af_threshold)
        & (hl.or_else(mt.gnomad_genomes.AF, MISSING_FLOAT_LO) < rare_af_threshold)
    )


def drop_useless_fields(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Remove fields from the MT which are irrelevant to the analysis
    Write times for the checkpoints are a substantial portion of the runtime
    Aim to reduce that by removing the amount of written data

    These fields would not be exported in the VCF anyway, so no downstream
    impacts caused by removal prior to that write

    Args:
        mt ():
    Returns:

    """

    # drop the useless top-level fields
    mt = mt.drop(*[field for field in USELESS_FIELDS if field in mt.row_value])

    # now drop most VEP fields
    return mt.annotate_rows(vep=hl.Struct(transcript_consequences=mt.vep.transcript_consequences))


def remove_variants_outside_gene_roi(mt: hl.MatrixTable, green_genes: hl.SetExpression) -> hl.MatrixTable:
    """
    chunky filtering method - get rid of every variant without at least one green-gene annotation
    does not edit the field itself, that will come later (split_rows_by_gene_and_filter_to_green)
    Args:
        mt ():
        green_genes ():

    Returns:
        the same MT, just reduced
    """

    # filter rows without a green gene (removes empty geneIds)
    return mt.filter_rows(hl.len(green_genes.intersection(mt.geneIds)) > 0)


def split_rows_by_gene_and_filter_to_green(mt: hl.MatrixTable, green_genes: hl.SetExpression) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any rows not annotating a Green PanelApp gene

    - first explode the matrix, separate gene per row
    - throw away all rows without a green gene
    - on all remaining rows, filter transcript consequences to match _this_ gene

    this is the single most powerful filtering row, effectively leaving just the genes we're interested in

    updated 25/10/2024 - we're retaining snRNA transcripts here, to enable them in

    Args:
        mt ():
        green_genes (): set of all relevant genes
    Returns:
        exploded MatrixTable
    """

    # split each gene onto a separate row, transforms 'geneIds' field from set to string
    mt = mt.explode_rows(mt.geneIds)

    # filter rows without a green gene (removes empty geneIds)
    mt = mt.filter_rows(green_genes.contains(mt.geneIds))

    # limit the per-row transcript CSQ to those relevant to the single
    # gene now present on each row
    return mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: (mt.geneIds == x.gene_id)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM') | (x.biotype == 'snRNA'))),
            ),
        ),
    )


def annotate_category_1(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Applies the boolean Category1 annotation
     - clinvar_talos_strong flag, as set in annotate_aip_clinvar
     - represents non-conflicting clinvar pathogenic/likely path

    Args:
        mt ():
    Returns:
        same variants, with categoryboolean1 set to 1 or 0
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean1=hl.if_else(mt.info.clinvar_talos_strong == ONE_INT, ONE_INT, MISSING_INT),
        ),
    )


def annotate_category_6(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the boolean Category6 flag
    - AlphaMissense likely Pathogenic on at least one transcript
    - Thresholds of am_pathogenicity:
        'Likely benign' if am_pathogenicity < 0.34;
        'Likely pathogenic' if am_pathogenicity > 0.564;
        'ambiguous' otherwise.

    If AM class attribute is missing, skip annotation (default to 0)
    This is run prior to Cat2 annotation - we use C6==True as a replacement
    for the previous approach of [CADD & REVEL predict damaging]

    Args:
        mt (hl.MatrixTable):
    Returns:
        same variants, categoryboolean6 set to 1 or 0
    """

    # focus on the auto-annotated AlphaMissense class
    # allow for the field to be missing
    if 'am_class' not in list(mt.vep.transcript_consequences[0].keys()):
        get_logger().warning('AlphaMissense class not found, skipping annotation')
        return mt.annotate_rows(info=mt.info.annotate(categoryboolean6=MISSING_INT))

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean6=hl.if_else(
                hl.len(mt.vep.transcript_consequences.filter(lambda x: x.am_class == 'likely_pathogenic')) > 0,
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def annotate_category_3(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the boolean Category3 flag
    - Critical protein consequence on at least one transcript
    - either predicted NMD or
    - any star Pathogenic or Likely_pathogenic in Clinvar

    Args:
        mt (hl.MatrixTable):
    Returns:
        same variants, categoryboolean3 set to 1 or 0
    """

    critical_consequences = hl.set(config_retrieve(['RunHailFiltering', 'critical_csq']))

    # First check if we have any HIGH consequences
    # then explicitly link the LOFTEE check with HIGH consequences
    # OR allow for a pathogenic ClinVar, any Stars
    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean3=hl.if_else(
                (
                    hl.len(
                        mt.vep.transcript_consequences.filter(
                            lambda x: (hl.len(critical_consequences.intersection(hl.set(x.consequence_terms))) > 0),
                        ),
                    )
                    > 0
                )
                & (
                    (
                        hl.len(
                            mt.vep.transcript_consequences.filter(
                                lambda x: (hl.len(critical_consequences.intersection(hl.set(x.consequence_terms))) > 0)
                                & ((x.lof == LOFTEE_HC) | (hl.is_missing(x.lof))),
                            ),
                        )
                        > 0
                    )
                    | (mt.info.clinvar_talos == ONE_INT)
                ),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def filter_by_consequence(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    - reduce the per-row transcript CSQ to a limited group
    - reduce the rows to ones where there are remaining tx consequences
    - alternatively except rows with spliceAI delta above threshold

    Args:
        mt ():
    """

    # at time of writing this is VEP HIGH + missense_variant
    # update without updating the dictionary content
    critical_consequences = set(config_retrieve(['RunHailFiltering', 'critical_csq'], [])) | set(
        config_retrieve(['RunHailFiltering', 'additional_csq'], []),
    )

    # overwrite the consequences with an intersection against a limited list
    filtered_mt = mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: (hl.len(hl.set(x.consequence_terms).intersection(critical_consequences)) > 0)
                | (x.biotype == 'snRNA'),
            ),
        ),
    )

    # filter out rows with no tx consequences left, and no splice cat. assignment
    return filtered_mt.filter_rows(
        (hl.len(filtered_mt.vep.transcript_consequences) == 0) & (filtered_mt.info.categoryboolean5 == 0),
        keep=False,
    )


def annotate_category_4(mt: hl.MatrixTable, ped_file_path: str) -> hl.MatrixTable:
    """
    Category based on de novo MOI, restricted to a group of consequences
    default uses the Hail builtin method (very strict)
    config switch to use the lenient version

    Args:
        mt ():
        ped_file_path (): path to a pedigree in PLINK format

    Returns:
        same variants, categorysample4 either 'missing' or sample IDs
        where de novo inheritance is seen
    """

    get_logger().info('Running de novo search')

    de_novo_matrix = filter_by_consequence(mt)

    pedigree = hl.Pedigree.read(ped_file_path)

    get_logger().info('Updating synthetic PL values for WT calls where missing')

    de_novo_matrix = de_novo_matrix.annotate_entries(
        PL=hl.case()
        .when(~hl.is_missing(de_novo_matrix.PL), de_novo_matrix.PL)
        .when((de_novo_matrix.GT.is_non_ref()) | (hl.is_missing(de_novo_matrix.GQ)), hl.missing('array<int32>'))
        .default([0, de_novo_matrix.GQ, 1000]),
    )

    dn_table = hl.de_novo(
        de_novo_matrix,
        pedigree,
        pop_frequency_prior=de_novo_matrix.info.gnomad_af,
        ignore_in_sample_allele_frequency=True,
        max_parent_ab=config_retrieve(['RunHailFiltering', 'max_parent_ab'], 0.05),
    )

    # re-key the table by locus,alleles, removing the sampleID from the compound key
    dn_table = dn_table.key_by(dn_table.locus, dn_table.alleles)

    # we only require the key (locus, alleles) and the sample ID
    # select to remove other fields, then collect per-key into Array of Structs
    dn_table = dn_table.select(dn_table.id).collect_by_key()

    # collect all sample IDs per locus, and squash into a String Array
    # delimit to compress that Array into single Strings
    dn_table = dn_table.annotate(dn_ids=hl.delimit(hl.map(lambda x: x.id, dn_table.values), ','))

    # log the number of variants found this way
    get_logger().info(f'{dn_table.count()} variants showed de novo inheritance')

    # annotate those values as a flag if relevant, else 'missing'
    return mt.annotate_rows(
        info=mt.info.annotate(categorysample4=hl.or_else(dn_table[mt.row_key].dn_ids, MISSING_STRING)),
    )


def annotate_category_5(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    SpliceAI based category assignment
    Args:
        mt ():

    Returns:
        same variants, categoryboolean5 set to 0 or 1
    """

    return mt.annotate_rows(
        info=mt.info.annotate(
            categoryboolean5=hl.if_else(
                mt.info.splice_ai_delta >= config_retrieve(['RunHailFiltering', 'spliceai']),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def vep_struct_to_csq(vep_expr: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a VEP Struct, returns an array of VEP VCF CSQ strings
    (1 per csq in the struct).
    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.
    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.

    Args:
        vep_expr (hl.Struct):
    Returns:
        generates an array of Strings for each CSQ
    """

    def get_csq_from_struct(element: hl.expr.StructExpression) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                'consequence': hl.delimit(element.consequence_terms, delimiter='&'),
                'feature': element.transcript_id,
                'ensp': element.protein_id,
                'gene': element.gene_id,
                'symbol': element.gene_symbol,
                'mane_select': element.mane_select,
            },
        )

        # pull the required fields and ordering from config
        csq_fields = config_retrieve(['RunHailFiltering', 'csq_string'])

        return hl.delimit([hl.or_else(hl.str(fields.get(f, '')), '') for f in csq_fields], '|')

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(vep_expr['transcript_consequences'].map(lambda x: get_csq_from_struct(x)), hl.empty_array(hl.tstr)),
    )

    # previous consequence filters may make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def filter_to_categorised(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter to rows tagged with a category label

    Args:
        mt ():
    Returns:
        input matrix, minus rows without Categories applied
    """

    return mt.filter_rows(
        (mt.info.categoryboolean1 == 1)
        | (mt.info.categoryboolean6 == 1)
        | (mt.info.categoryboolean3 == 1)
        | (mt.info.categorysample4 != MISSING_STRING)
        | (mt.info.categoryboolean5 == 1)
        | (mt.info.categorydetailspm5 != MISSING_STRING)
        | (mt.info.categorybooleansvdb == 1)
        | (mt.info.categorydetailsexomiser != MISSING_STRING),
    )


def write_matrix_to_vcf(mt: hl.MatrixTable, vcf_out: str):
    """
    write the remaining MatrixTable content to file as a VCF
    generate a custom header containing the CSQ contents which
    were retained during this run

    Args:
        mt (): the whole MatrixTable
        vcf_out (str): where to write the VCF
    """

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    header_path = 'additional_header.txt'

    # generate a CSQ string specific to the config file for decoding later
    csq_contents = '|'.join(config_retrieve(['RunHailFiltering', 'csq_string']))

    # write this custom header locally
    with open(header_path, 'w') as handle:
        handle.write(f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {csq_contents}">')
    get_logger().info(f'Writing categorised variants out to {vcf_out}')
    hl.export_vcf(mt, vcf_out, append_to_header=header_path, tabix=True)


def green_from_panelapp(panel_data: PanelApp) -> hl.SetExpression:
    """
    Pull all ENSGs from PanelApp data relating to Green Genes

    Args:
        panel_data (PanelApp): the PanelApp object

    Returns:
        a set expression - all Green genes in the analysis
    """

    # take all the green genes, remove the metadata
    green_genes = set(panel_data.genes.keys())
    get_logger().info(f'Extracted {len(green_genes)} green genes')
    return hl.literal(green_genes)


def subselect_mt_to_pedigree(mt: hl.MatrixTable, pedigree: str) -> hl.MatrixTable:
    """
    remove any columns from the MT which are not represented in the Pedigree
    _not_ done: remove variants without calls amongst the current samples
    that's probably more processing than benefit - will be skipped later

    Args:
        mt ():
        pedigree ():
    Returns:
        MatrixTable, possibly with fewer columns (samples) present
    """

    # individual IDs from pedigree
    ped_samples: set[str] = set()
    for family in open_ped(pedigree):
        ped_samples.update({member.id for member in family})

    # individual IDs from matrix
    matrix_samples = set(mt.s.collect())

    # find overlapping samples
    common_samples: set[str] = ped_samples.intersection(matrix_samples)

    get_logger().info(
        f"""
    Samples in Pedigree: {len(ped_samples)}
    Samples in MatrixTable: {len(matrix_samples)}
    Common Samples: {len(common_samples)}
    """,
    )

    if not common_samples:
        raise ValueError('No samples shared between pedigree and MT')

    # full overlap = no filtering
    if common_samples == matrix_samples:
        return mt

    # reduce to those common samples
    mt = mt.filter_cols(hl.literal(common_samples).contains(mt.s))
    get_logger().info(f'Remaining MatrixTable columns: {mt.count_cols()}')

    return mt


def generate_a_checkpoint(mt: hl.MatrixTable, checkpoint_path: str) -> hl.MatrixTable:
    """
    wrapper around a few repeated lines which write a checkpoint
    Args:
        mt ():
        checkpoint_path (str): where to write the checkpoint

    Returns:

    """
    get_logger().info(f'Checkpointing to {checkpoint_path} after filtering out a ton of variants')
    mt = mt.checkpoint(checkpoint_path)

    # die if there are no variants remaining. Only ever count rows after a checkpoint
    if not (current_rows := mt.count_rows()):
        raise ValueError('No remaining rows to process!')

    get_logger().info(f'Local checkpoint written, {current_rows} rows remain')
    return mt


def cli_main():
    """
    Read MT, filter, and apply category annotation, export as a VCF
    """

    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='path to input MT')
    parser.add_argument('--panelapp', required=True, help='panelapp JSON')
    parser.add_argument('--pedigree', required=True, help='Cohort Pedigree')
    parser.add_argument('--output', help='Where to write the VCF', required=True)
    parser.add_argument('--clinvar', help='HT containing ClinvArbitration annotations', required=True)
    parser.add_argument('--pm5', help='HT containing clinvar PM5 annotations, optional', default=None)
    parser.add_argument('--exomiser', help='HT containing exomiser variant selections, optional', default=None)
    parser.add_argument('--svdb', help='HT containing SpliceVarDB annotations, optional', default=None)
    parser.add_argument('--checkpoint', help='Where/whether to checkpoint, String path', default=None)
    args = parser.parse_args()
    main(
        mt_path=args.input,
        panel_data=args.panelapp,
        pedigree=args.pedigree,
        vcf_out=args.output,
        clinvar=args.clinvar,
        pm5=args.pm5,
        exomiser=args.exomiser,
        svdb=args.svdb,
        checkpoint=args.checkpoint,
    )


def main(
    mt_path: str,
    panel_data: str,
    pedigree: str,
    vcf_out: str,
    clinvar: str,
    pm5: str | None = None,
    exomiser: str | None = None,
    svdb: str | None = None,
    checkpoint: str | None = None,
):
    """
    Read MT, filter, and apply category annotation, export as a VCF

    Args:
        mt_path (str): Location of the input MT
        panel_data ():
        pedigree (str): path to a pedigree for this cohort (used in de novo testing)
        vcf_out (str): where to write VCF out
        clinvar (str): path to a ClinVar HT, or unspecified
        pm5 (str): path to a pm5 HT, or unspecified
        exomiser (str): path of an exomiser HT, or unspecified
        svdb (str): path to a SpliceVarDB HT, or unspecified
        checkpoint (str): path to checkpoint data to - serves as checkpoint trigger
    """
    get_logger(__file__).info(
        r"""Welcome To
███████████   █████████   █████          ███████     █████████
█   ███   █  ███     ███   ███         ███     ███  ███     ███
    ███      ███     ███   ███        ███       ███ ███
    ███      ███████████   ███        ███       ███  █████████
    ███      ███     ███   ███        ███       ███         ███
    ███      ███     ███   ███      █  ███     ███  ███     ███
   █████    █████   █████ ███████████    ███████     █████████""",
    )

    # initiate Hail as a local cluster
    number_of_cores = config_retrieve(['RunHailFiltering', 'cores', 'small_variants'], 8)
    get_logger().info(f'Starting Hail with reference genome GRCh38, as a {number_of_cores} core local cluster')

    hl.context.init_spark(master=f'local[{number_of_cores}]', default_reference='GRCh38', quiet=True)

    # read the parsed panelapp data
    get_logger().info(f'Reading PanelApp data from {panel_data!r}')
    panelapp = read_json_from_path(panel_data, return_model=PanelApp)
    assert isinstance(panelapp, PanelApp)

    # pull green genes from the panelapp data
    green_expression = green_from_panelapp(panelapp)

    # read the matrix table from a localised directory
    mt = hl.read_matrix_table(mt_path)
    get_logger().info(f'Loaded annotated MT from {mt_path}, size: {mt.count_rows()}, partitions: {mt.n_partitions()}')

    # lookups for required fields all delegated to the hail_audit file
    if not (
        fields_audit(mt=mt, base_fields=BASE_FIELDS_REQUIRED, nested_fields=FIELDS_REQUIRED)
        and vep_audit(mt=mt, expected_fields=VEP_TX_FIELDS_REQUIRED)
    ):
        raise KeyError('Fields were missing from the input Matrix')

    # repartition if required - local Hail with finite resources has struggled with some really high (~120k) partitions
    # this creates a local duplicate of the input data with far smaller partition counts, for less processing overhead
    if mt.n_partitions() > MAX_PARTITIONS:
        get_logger().info('Shrinking partitions way down with a unshuffled repartition')
        mt = mt.repartition(shuffle=False, n_partitions=number_of_cores * 10)
        if checkpoint:
            get_logger().info('Trying to write the result locally, might need more space on disk...')
            mt = generate_a_checkpoint(mt, f'{checkpoint}_reparitioned')

    # subset to currently considered samples
    mt = subselect_mt_to_pedigree(mt, pedigree=pedigree)

    # shrink the time taken to write checkpoints
    mt = drop_useless_fields(mt=mt)

    # remove any rows which have no genes of interest
    mt = remove_variants_outside_gene_roi(mt=mt, green_genes=green_expression)

    # remove common-in-gnomad variants
    mt = filter_to_population_rare(mt=mt)

    if checkpoint:
        mt = generate_a_checkpoint(mt, f'{checkpoint}_green_genes')

    # swap out the default clinvar annotations with private clinvar
    mt = annotate_clinvarbitration(mt=mt, clinvar=clinvar)

    # annotate this MT with exomiser variants - annotated as MISSING if the table is absent
    mt = annotate_exomiser(mt=mt, exomiser=exomiser)

    # if a SVDB data is provided, use that to apply category annotations
    mt = annotate_splicevardb(mt=mt, svdb_path=svdb)

    # filter out quality failures
    mt = filter_on_quality_flags(mt=mt)

    # running global quality filter steps
    mt = filter_to_well_normalised(mt=mt)

    # filter variants by frequency - allow for ClinVar as an escape mechanism for this filter if
    if config_retrieve(['RunHailFiltering', 'allow_clinvar_common'], False):
        mt = permissive_filter_matrix_by_ac(mt=mt)
    else:
        mt = filter_matrix_by_ac(mt=mt)

    # rearrange the row annotation to make syntax nicer downstream
    mt = extract_annotations(mt=mt)

    # split each gene annotation onto separate rows, filter to green genes (PanelApp ROI)
    mt = split_rows_by_gene_and_filter_to_green(mt=mt, green_genes=green_expression)

    if checkpoint:
        mt = generate_a_checkpoint(mt, f'{checkpoint}_green_and_clean')

    # add Labels to the MT
    # current logic is to apply 1, 2, 3, and 5, then 4 (de novo)
    # for cat. 4, pre-filter the variants by tx-consequential or C5==1
    get_logger().info('Applying categories')
    mt = annotate_category_1(mt=mt)
    mt = annotate_category_6(mt=mt)
    mt = annotate_category_3(mt=mt)
    mt = annotate_category_5(mt=mt)

    # ordering is important - category4 (de novo) makes
    # use of category 5, so it must follow
    mt = annotate_category_4(mt=mt, ped_file_path=pedigree)

    # if a clinvar-codon table is supplied, use that for PM5
    mt = annotate_codon_clinvar(mt=mt, pm5_path=pm5)

    # remove all variants without positively assigned labels
    mt = filter_to_categorised(mt=mt)

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    mt = mt.annotate_rows(info=mt.info.annotate(CSQ=vep_struct_to_csq(mt.vep), gene_id=mt.geneIds))

    write_matrix_to_vcf(mt=mt, vcf_out=vcf_out)


if __name__ == '__main__':
    cli_main()
