
process ValidateMOI {
    container params.container


    // process the labelled variants
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple path(labelled_vcf), path(labelled_vcf_index)
        path panelapp
        path pedigree
        path hpo_panel_matches
        path talos_config

    output:
        path"${params.cohort}_results.json"

    """
    export TALOS_CONFIG=${talos_config}

    ValidateMOI \
        --labelled_vcf ${labelled_vcf} \
        --panelapp ${panelapp} \
        --pedigree ${pedigree} \
        --participant_panels ${hpo_panel_matches} \
        --output ${params.cohort}_results.json
    """
}
