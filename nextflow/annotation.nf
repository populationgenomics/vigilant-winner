#!/usr/bin/env nextflow

/*
This workflow is a minimised annotation process for the Talos pipeline.

It takes multiple VCFs, merges them into a single VCF, and annotates the merged VCF with relevant annotations.

The specific annotations are:
- gnomAD v4.1 frequencies, applied to the joint VCF using echtvar
- Transcript consequences, using BCFtools annotate
- AlphaMissense, applied using Hail
- MANE trancript IDs and corresponding proteins, applied using Hail
*/

nextflow.enable.dsl=2

include { ParseAlphaMissenseIntoHt } from './modules/annotation/ParseAlphaMissenseIntoHt/main'
include { AnnotateCsqWithBcftools } from './modules/annotation/AnnotateCsqWithBcftools/main'
include { CreateRoiFromGff3 } from './modules/annotation/CreateRoiFromGff3/main'
include { AnnotateGnomadAfWithEchtvar } from './modules/annotation/AnnotateGnomadAfWithEchtvar/main'
include { LocaliseAlphamissenseWithWget } from './modules/annotation/LocaliseAlphamissenseWithWget/main'
include { MergeVcfsWithBcftools } from './modules/annotation/MergeVcfsWithBcftools/main'
include { ParseManeIntoJson } from './modules/annotation/ParseManeIntoJson/main'
include { ReformatVcfToMt } from './modules/annotation/ReformatVcfToMt/main'

workflow {
    // generate the AlphaMissense HT - long running, stored in a separate folder
    if (file(params.alphamissense_output).exists()) {
        ch_alphamissense_table = channel.fromPath(params.alphamissense_output)
    }
    else {
        LocaliseAlphamissenseWithWget()
        ParseAlphaMissenseIntoHt(LocaliseAlphamissenseWithWget.out)
        ch_alphamissense_table = ParseAlphaMissenseIntoHt.out
    }

    // generate the gene region file, and a overlap-merged version of the same
    CreateRoiFromGff3()

    // get all the VCFs
    ch_vcfs = channel.fromPath(params.input_vcfs)
    ch_tbis = channel.fromPath(params.input_vcfs).map{ it -> file("${it}.tbi") }

    MergeVcfsWithBcftools(
        ch_vcfs.collect(),
        ch_tbis.collect(),
        CreateRoiFromGff3.out.merged_bed,
    )

    // read the whole-genome Zip file as an input channel
    ch_gnomad_zip = channel.fromPath(params.gnomad_zip)

    // and apply the annotation using echtvar
    AnnotateGnomadAfWithEchtvar(
        MergeVcfsWithBcftools.out,
        ch_gnomad_zip,
    )

    // bcftools csq
    ch_ref_genome = channel.fromPath(params.ref_genome)
    AnnotateCsqWithBcftools(
        AnnotateGnomadAfWithEchtvar.out,
        CreateRoiFromGff3.out.gff3,
        ch_ref_genome
    )

    // pull and parse the MANE data into a Hail Table
    ParseManeIntoJson()

    // now what
    ReformatVcfToMt(
        AnnotateCsqWithBcftools.out,
        ch_alphamissense_table,
        CreateRoiFromGff3.out.bed,
        ParseManeIntoJson.out.json
    )
}
