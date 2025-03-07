
process AnnotateGnomadAfWithEchtvar {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path(zip)

    // annotate this VCF with gnomAD data
    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}_gnomad.vcf.bgz")

    script:
        """
        set -ex
        sorted_zips=""

        for chrom in \$(seq 1 22) X Y; do
            sorted_zips="\${sorted_zips} -e chr\${chrom}.zip "
        done;

        echtvar anno \
            -e ${zip} \
            ${vcf} \
            "${params.cohort}_gnomad.vcf.bgz"
        """
}
