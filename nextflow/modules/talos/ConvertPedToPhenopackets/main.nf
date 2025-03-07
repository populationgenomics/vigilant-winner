
process ConvertPedToPhenopackets {
    container params.container

    // takes the pedigree file, and converts it to a phenopackets file and regular pedigree file
    publishDir params.output_dir, mode: 'copy'

    input:
        // the pedigree with embedded HPO terms
        path pedigree

    output:
        path "${params.cohort}_pedigree.ped", emit: "ped"
        path "${params.cohort}_phenopackets.json", emit: "phenopackets"

    """
    ConvertPedToPhenopackets \
        --input ${pedigree} \
        --output ${params.cohort}
    """
}
