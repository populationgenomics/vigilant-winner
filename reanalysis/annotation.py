"""
Hail Query function to convert VCF to a MatrixTable, and add annotations
based on reference data: VEP, ClinVar, etc.
"""


import logging

from cloudpathlib import AnyPath
from cpg_utils.hail_batch import reference_path
import hail as hl


logger = logging.getLogger(__file__)


def apply_annotations(
    vcf_path: str,
    vep_ht_path: str,
    out_mt_path: str,
    overwrite: bool = False,
    genome_build: str = 'GRCh38',
    sequencing_type: str = 'WGS',
    checkpoints_bucket=None,
    ref_ht_path=None,
    clinvar_ht_path=None,
):
    """
    Hail Query function to convert VCF to a MatrixTable, and add annotations
    based on reference data: VEP, ClinVar, etc.
    """
    # pylint: disable=C0415
    from hail_scripts.computed_fields.variant_id import (
        get_expr_for_contig,
        get_expr_for_variant_id,
        get_expr_for_variant_ids,
        get_expr_for_xpos,
    )

    # custom imports to reduce line lengths
    # pylint disable:C0415
    from hail_scripts.computed_fields.vep import (
        get_expr_for_vep_sorted_transcript_consequences_array as vep_array,
        get_expr_for_vep_transcript_ids_set as vep_set,
        get_expr_for_worst_transcript_consequence_annotations_struct as worst_csq,
        get_expr_for_vep_consequence_terms_set as vep_csq_set,
        get_expr_for_vep_gene_ids_set as vep_gene_id_set,
        get_expr_for_vep_protein_domains_set_from_sorted as vep_prot,
    )

    ref_ht_path = ref_ht_path or reference_path('seqr/combined_reference')
    clinvar_ht_path = clinvar_ht_path or reference_path('seqr/clinvar')

    assert AnyPath(clinvar_ht_path).exists(), f'{clinvar_ht_path} unavailable'
    assert AnyPath(ref_ht_path).exists(), f'{ref_ht_path} unavailable'

    def _checkpoint(t: hl.Table, filename: str):
        """
        if a checkpoint location was provided, periodically dump the table to
        that location

        :param t:
        :param filename:
        :return:
        """
        if checkpoints_bucket:
            t = t.checkpoint(
                f'{checkpoints_bucket}/{filename}', _read_if_exists=not overwrite
            )
            logger.info(f'Wrote {checkpoints_bucket}/{filename}')
        return t

    mt = hl.import_vcf(
        str(vcf_path),
        reference_genome=genome_build,
        skip_invalid_loci=True,
        force_bgz=True,
    )
    logger.info(
        f'Importing VCF {vcf_path}, ' f'adding VEP annotations from {vep_ht_path}'
    )

    logger.info(f'Loading VEP Table from {vep_ht_path}')
    # Annotate VEP. Do before splitting multi, because we run VEP on unsplit VCF,
    # and hl.split_multi_hts can handle multiallelic VEP field.
    vep_ht = hl.read_table(str(vep_ht_path))
    logger.info(f'Adding VEP annotations into the Matrix Table from {vep_ht_path}')
    mt = mt.annotate_rows(vep=vep_ht[mt.locus].vep)

    # Splitting multi-allelics. We do not handle AS info fields here
    mt = hl.split_multi_hts(
        mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles)
    )
    mt = _checkpoint(mt, 'mt-vep-split.mt')

    # Add dummy GRCh37 loci
    logger.info('Adding dummy GRCh37 loci')
    mt = mt.annotate_rows(rg37_locus=mt.locus)

    ref_ht = hl.read_table(str(ref_ht_path))
    clinvar_ht = hl.read_table(str(clinvar_ht_path))

    logger.info('Annotating with seqr-loader fields: round 1')
    logger.info(mt.describe())
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(
        info=mt.info.annotate(AN=mt.variant_qc.AN, AF=mt.variant_qc.AF)
    )
    mt = mt.drop('variant_qc')
    mt = mt.annotate_rows(
        AC=mt.info.AC,
        AF=mt.info.AF[mt.a_index - 1],
        AN=mt.info.AN,
        aIndex=mt.a_index,
        wasSplit=mt.was_split,
        originalAltAlleles=get_expr_for_variant_ids(mt.locus_old, mt.alleles_old),
        sortedTranscriptConsequences=vep_array(mt.vep),
        variantId=get_expr_for_variant_id(mt),
        contig=get_expr_for_contig(mt.locus),
        pos=mt.locus.position,
        start=mt.locus.position,
        end=mt.locus.position + hl.len(mt.alleles[0]) - 1,
        ref=mt.alleles[0],
        alt=mt.alleles[1],
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(mt.locus) + hl.len(mt.alleles[0]) - 1,
        rg37_locus=mt.rg37_locus,
        clinvar_data=clinvar_ht[mt.row_key],
        ref_data=ref_ht[mt.row_key],
    )
    mt = _checkpoint(mt, 'mt-vep-split-round1.mt')

    logger.info(
        'Annotating with seqr-loader fields: round 2 '
        '(expanding sortedTranscriptConsequences, ref_data, clinvar_data)'
    )
    mt = mt.annotate_rows(
        domains=vep_prot(mt.sortedTranscriptConsequences),
        transcriptConsequenceTerms=vep_csq_set(mt.sortedTranscriptConsequences),
        transcriptIds=vep_set(mt.sortedTranscriptConsequences),
        mainTranscript=worst_csq(mt.sortedTranscriptConsequences),
        geneIds=vep_gene_id_set(mt.sortedTranscriptConsequences),
        codingGeneIds=vep_gene_id_set(
            mt.sortedTranscriptConsequences, only_coding_genes=True
        ),
        cadd=mt.ref_data.cadd,
        dbnsfp=mt.ref_data.dbnsfp,
        geno2mp=mt.ref_data.geno2mp,
        gnomad_exomes=mt.ref_data.gnomad_exomes,
        gnomad_exome_coverage=mt.ref_data.gnomad_exome_coverage,
        gnomad_genomes=mt.ref_data.gnomad_genomes,
        gnomad_genome_coverage=mt.ref_data.gnomad_genome_coverage,
        eigen=mt.ref_data.eigen,
        exac=mt.ref_data.exac,
        g1k=mt.ref_data.g1k,
        mpc=mt.ref_data.mpc,
        primate_ai=mt.ref_data.primate_ai,
        splice_ai=mt.ref_data.splice_ai,
        topmed=mt.ref_data.topmed,
        clinvar=hl.struct(
            **{
                'allele_id': mt.clinvar_data.info.ALLELEID,
                'clinical_significance': hl.delimit(mt.clinvar_data.info.CLNSIG),
                'gold_stars': mt.clinvar_data.gold_stars,
            }
        ),
    )
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=genome_build.replace('GRCh', ''),
        sampleType=sequencing_type,
        hail_version=hl.version(),
    )
    mt.write(str(out_mt_path), overwrite=overwrite)
    logger.info(f'Written final matrix table into {out_mt_path}')
