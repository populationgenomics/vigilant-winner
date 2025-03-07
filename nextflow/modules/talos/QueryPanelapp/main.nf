
process QueryPanelapp {
    container params.container


    // uses matched panels to query the PanelApp API
    publishDir params.output_dir, mode: 'copy'

    input:
        path hpo_panel_matches
        path talos_config

    output:
        path "${params.cohort}_panelapp_results.json"

    // the command
    """
    export TALOS_CONFIG=${talos_config}
    QueryPanelapp --input ${hpo_panel_matches}  --output ${params.cohort}_panelapp_results.json
    """
}
