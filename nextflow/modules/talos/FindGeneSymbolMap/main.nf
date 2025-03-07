
process FindGeneSymbolMap {
    container params.container

    // finds corresponding gene symbols for the panelapp results
    publishDir params.output_dir, mode: 'copy'

    input:
        path panelapp_data
        path talos_config

    output:
        path "${params.cohort}_symbol_to_ensg.json"

    // the command
    """
    export TALOS_CONFIG=${talos_config}
    FindGeneSymbolMap --input ${panelapp_data}  --output ${params.cohort}_symbol_to_ensg.json
    """
}
