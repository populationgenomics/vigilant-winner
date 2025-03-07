
process ConvertSpliceVarDb {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        // input splicevardb TSV
        path svdb

    output:
        path "${params.cohort}_svdb.ht.tar.gz"

    """
    ConvertSpliceVarDb \
        --input ${svdb} \
        --output ${params.cohort}_svdb.ht
    tar -czf ${params.cohort}_svdb.ht.tar.gz ${params.cohort}_svdb.ht
    """
}
