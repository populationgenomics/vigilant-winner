
process ParseAlphaMissenseIntoHt {
    container params.container

    // parse AM data as a Hail Table
    publishDir params.generic_output_dir, mode: 'copy'

    input:
        path(am_tsv)

    output:
        path("alphamissense_38.ht.tar.gz")

    script:
        """
        ParseAlphaMissenseIntoHt \
            --am_tsv ${am_tsv} \
            --ht_out alphamissense_38.ht
        tar -czf alphamissense_38.ht.tar.gz alphamissense_38.ht
        rm -r alphamissense_38.ht
        """
}
