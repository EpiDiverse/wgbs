#!/usr/bin/env nextflow


// FUNCTION TO TEST REFERENCE PARAMS IN EPI PROFILE
def check_ref_errors(reference,thlaspi,fragaria,populus,nolambda) {

    // PATH TO BASE DIRECTORIES OF EPIDIVERSE SPECIES
    def thlaspi_dir = "/scr/epi/genomes/thlaspi_arvense"
    def fragari_dir = "/scr/epi/genomes/fragaria_vesca"
    def populus_dir = "/scr/epi/genomes/populus_nigra"

    // error if multiple references given
    if (((( thlaspi ) && (( fragaria ) || ( populus ) || ( reference ))) || \
    (( fragaria ) && (( populus ) || ( reference ))) || (( populus ) && ( reference )))) {
        error "ERROR: more than one reference genome has been specified!"
        exit 1
    }

    // params.reference
    else if ( reference ) {
        def fasta = file("${params.reference}", checkIfExists: true, glob: false)
        def fai = file("${params.reference}.fai", checkIfExists: true, glob: false)
        def ebm = nolambda ? "${fasta.parent}/index/*.ebm" : "${fasta.parent}/lambda/*.ebm" 
        def ctidx = nolambda ? "${fasta.parent}/index/*.ctidx" : "${fasta.parent}/lambda/*.ctidx"
        def gaidx = nolambda ? "${fasta.parent}/index/*.gaidx" : "${fasta.parent}/lambda/*.gaidx"
        return tuple(fasta, fai, ebm, ctidx, gaidx)
    }

    // params.thlaspi
    else if ( params.thlaspi ) {
        def fasta = file("${thlaspi_dir}/thlaspi.fa", checkIfExists: true, glob: false)
        def fai = file("${thlaspi_dir}/thlaspi.fa.fai", checkIfExists: true, glob: false)
        def ebm = nolambda ? "${thlaspi_dir}/index/thlaspi.ebm" : "${thlaspi_dir}/lambda/lambda.ebm"
        def ctidx = nolambda ? "${thlaspi_dir}/index/thlaspi.ctidx" : "${thlaspi_dir}/lambda/lambda.ctidx"
        def gaidx = nolambda ? "${thlaspi_dir}/index/thlaspi.gaidx" : "${thlaspi_dir}/lambda/lambda.gaidx"
        return tuple(fasta, fai, ebm, ctidx, gaidx)
    }

    // params.fragaria
    else if ( fragaria ) {
        def fasta = file("${fragari_dir}/fragaria.fa", checkIfExists: true, glob: false)
        def fai = file("${fragari_dir}/fragaria.fa.fai", checkIfExists: true, glob: false)
        def ebm = nolambda ? "${fragari_dir}/index/fragaria.ebm" : "${fragari_dir}/lambda/lambda.ebm"
        def ctidx = nolambda ? "${fragari_dir}/index/fragaria.ctidx" : "${fragari_dir}/lambda/lambda.ctidx"
        def gaidx = nolambda ? "${fragari_dir}/index/fragaria.gaidx" : "${fragari_dir}/lambda/lambda.gaidx"
        return tuple(fasta, fai, ebm, ctidx, gaidx)
    }

    // params.populus
    else if ( populus ) {
        def fasta = file("${populus_dir}/populus.fa", checkIfExists: true, glob: false)
        def fai = file("${populus_dir}/populus.fa.fai", checkIfExists: true, glob: false)
        def ebm = nolambda ? "${populus_dir}/index/populus.ebm" : "${populus_dir}/lambda/lambda.ebm"
        def ctidx = nolambda ? "${populus_dir}/index/populus.ctidx" : "${populus_dir}/lambda/lambda.ctidx"
        def gaidx = nolambda ? "${populus_dir}/index/populus.gaidx" : "${populus_dir}/lambda/lambda.gaidx"
        return tuple(fasta, ebm, ctidx, gaidx)
    }

    else {
        error "ERROR: please specify a reference genome."
        exit 1
    }
}


// FUNCTION TO LOAD DATASETS IN TEST PROFILE
def check_test_data(readPaths,mergePaths,singleEnd,merge) {

    // SINGLE END TESTDATA
    if( singleEnd ){
        reads = Channel
            .from(readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        merged = Channel
            .from(mergePaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.mergePaths was empty - no input files supplied" }

    // PAIRED END TESTDATA
    } else {

        reads = Channel
            .from(readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        merged = Channel
            .from(mergePaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.mergePaths was empty - no input files supplied" }
    }

    // return the reads from testdata repo
    if( merge ) {
        return tuple(reads, merged)
    } else {
        return tuple(reads, Channel.empty())
    }

}