

process AnnotateCsqWithBcftools {
    container params.container

    input:
        path(vcf)
        path(gff3)
        path(reference)

    // annotate this VCF with gnomAD data
    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}_csq.vcf.bgz")

    script:
    """
    bcftools index -t ${vcf}
    bcftools csq --force -f "${reference}" \
        --local-csq \
        -g ${gff3} \
        -Oz -o "${params.cohort}_csq.vcf.bgz" \
        ${vcf}
    """
}
