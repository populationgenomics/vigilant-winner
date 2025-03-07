
process GeneratePanelData {
    container params.container


    // takes the HPO-embellished pedigree and matches panels to participants
    publishDir params.output_dir, mode: 'copy'

    input:
        // the pedigree file
        path pedigree
        // the HPO obo/ontology file
        path hpo

    output:
        path "${params.cohort}_hpo_panel_data.json"

    """
    GeneratePanelData --input ${pedigree} --hpo ${hpo} --output ${params.cohort}_hpo_panel_data.json
    """
}
