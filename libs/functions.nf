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
        def fasta = "${params.reference}/*.fa"
        def ebm = nolambda ? "${reference}/index/*.ebm" : "${reference}/lambda/*.ebm" 
        def ctidx = nolambda ? "${reference}/index/*.ctidx" : "${reference}/lambda/*.ctidx"
        def gaidx = nolambda ? "${reference}/index/*.gaidx" : "${reference}/lambda/*.gaidx"
        return tuple(fasta, ebm, ctidx, gaidx)
    }

    // params.thlaspi
    else if ( params.thlaspi ) {
        def fasta = "${thlaspi_dir}/thlaspi.fa" 
        def ebm = nolambda ? "${thlaspi_dir}/index/thlaspi.ebm" : "${thlaspi_dir}/lambda/lambda.ebm"
        def ctidx = nolambda ? "${thlaspi_dir}/index/thlaspi.ctidx" : "${thlaspi_dir}/lambda/lambda.ctidx"
        def gaidx = nolambda ? "${thlaspi_dir}/index/thlaspi.gaidx" : "${thlaspi_dir}/lambda/lambda.gaidx"
        return tuple(fasta, ebm, ctidx, gaidx)
    }

    // params.fragaria
    else if ( fragaria ) {
        def fasta = "${fragari_dir}/fragaria.fa"
        def ebm = nolambda ? "${fragari_dir}/index/fragaria.ebm" : "${fragari_dir}/lambda/lambda.ebm"
        def ctidx = nolambda ? "${fragari_dir}/index/fragaria.ctidx" : "${fragari_dir}/lambda/lambda.ctidx"
        def gaidx = nolambda ? "${fragari_dir}/index/fragaria.gaidx" : "${fragari_dir}/lambda/lambda.gaidx"
        return tuple(fasta, ebm, ctidx, gaidx)
    }

    // params.populus
    else if ( populus ) {
        def fasta = "${populus_dir}/populus.fa"
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