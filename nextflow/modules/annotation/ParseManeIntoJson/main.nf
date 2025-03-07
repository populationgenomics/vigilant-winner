
process ParseManeIntoJson {
    container params.container

    publishDir params.generic_output_dir, mode: 'copy'

    output:
        path "mane_summary.txt.gz", emit: summary
        path "mane.json", emit: json

    script:
    """
    wget ${params.mane} -O mane_summary.txt.gz
    ParseManeIntoJson --input mane_summary.txt.gz --output mane.json --format json
    """
}
