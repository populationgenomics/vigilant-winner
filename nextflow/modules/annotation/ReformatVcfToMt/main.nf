
process ReformatVcfToMt {
    container params.container

    input:
        path(vcf)
        path(alphamissense)
        path(gene_bed)
        path(mane)

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}.mt.tar.gz")

    script:
        """
        set -ex
        tar -xf ${alphamissense}
        ReformatVcfToMt \
            --input ${vcf} \
            --am alphamissense_38.ht \
            --gene_bed ${gene_bed} \
            --output ${params.cohort}.mt \
            --mane ${mane}
        # cut down on work folder space
        rm -r alphamissense_38.ht

        tar -czf ${params.cohort}.mt.tar.gz ${params.cohort}.mt
        # not deleting the MT at the moment, someone might want to use it without decompressing the process output
        """
}
