
process LocaliseAlphamissenseWithWget {
    container params.container

    // parse AM data as a Hail Table
    publishDir params.generic_output_dir, mode: 'copy'

    output:
        path("alphamissense_38.tsv.gz")

    script:
        """
        wget ${params.alphamissense_url} -O alphamissense_38.tsv.gz
        """
}
