
process RunHailFiltering {
    container params.container

    // runs the hail small-variant filtering
    publishDir params.output_dir, mode: 'copy'

    input:
        path mt_tar
        path panelapp_data
        path pedigree
        path clinvar
        path exomiser
        path svdb
        val checkpoint
        path talos_config

    output:
        tuple \
            path("${params.cohort}_small_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_small_variants_labelled.vcf.bgz.tbi")

    // only write a checkpoint if we were given a path
    script:
        def checkpoint = checkpoint != "${params.output_dir}/assets/NO_FILE" ? "--checkpoint ${checkpoint}" : ''

    // unzip the ClinvArbitration data directory and MatrixTable
    """
    export TALOS_CONFIG=${talos_config}

    tar -zxf ${clinvar}
    tar -zxf ${mt_tar}
    tar -zxf ${svdb}
    tar -zxf ${exomiser}

    RunHailFiltering \
        --input ${params.cohort}.mt \
        --panelapp ${panelapp_data} \
        --pedigree ${pedigree} \
        --output ${params.cohort}_small_variants_labelled.vcf.bgz \
        --clinvar clinvarbitration_data/clinvar_decisions.ht \
        --pm5 clinvarbitration_data/clinvar_pm5.ht \
        --svdb ${params.cohort}_svdb.ht \
        --exomiser exomiser_nf.ht \
        ${checkpoint}
    """
}
