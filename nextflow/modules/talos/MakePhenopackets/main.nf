
process MakePhenopackets {
    container params.container


    publishDir params.output_dir, mode: 'copy'

    input:
        // the name of the dataset for use in API queries
        val dataset

        // the sequencing type
        val sequencing_type

        // the path to the HPO file
        path hpo

        // sequencing technology, not currently overridden
        val sequencing_tech

    output:
        path "${params.cohort}_phenopackets.json"

    """
    MakePhenopackets \
        --dataset ${dataset} \
        --output ${params.cohort}_phenopackets.json \
        --type ${sequencing_type} \
        --hpo ${hpo} \
        --tech ${sequencing_tech}
    """
}
