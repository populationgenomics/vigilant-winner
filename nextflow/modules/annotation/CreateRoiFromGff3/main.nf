
process CreateRoiFromGff3 {
    container params.container

    // generate the ROI file
    publishDir params.generic_output_dir, mode: 'copy'

    // generates a BED file one-gene-per-line, for per-gene annotations
    // generates an BED with overlapping intervals merged for filtering
    // sorts both
    // retains the downloaded GFF3 file
    output:
        path "GRCh38.bed", emit: bed
        path "GRCh38_merged.bed", emit: merged_bed
        path "GRCh38.gff3.gz", emit: gff3

    script:
        def url_template = "${params.ensembl_gtf}".replaceAll("VER", "${params.ensembl_version}")
        """
        wget ${url_template} -O GRCh38.gff3.gz
        CreateRoiFromGff3 \
            --gff3 GRCh38.gff3.gz \
            --unmerged_output unsorted_GRCh38.bed \
            --merged_output unsorted_merged_GRCh38.bed

        sort -k1,1V -k2,2n unsorted_GRCh38.bed > GRCh38.bed
        sort -k1,1V -k2,2n unsorted_merged_GRCh38.bed > GRCh38_merged.bed
        """
}
