#!/usr/bin/env nextflow

// DSL2 BRANCH - this branch is for testing new Nextflow features
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         ===============================================
         E P I D I V E R S E - W G B S   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run epidiverse/wgbs [OPTIONS]...


         Options: INPUT/OUTPUT
              --input [path/to/reads/dir]     [REQUIRED] Specify the path to the DIRECTORY containing the reads for analysis.
                                          All fastq files must be in compressed '*{1,2}.fastq.gz' format, with separate files
                                          for R1 and R2. This format must be used for both PE and SE reads, but these cannot
                                          be processed together. Any number of read sets (of the same type) may be included,
                                          and they will be analysed in parallel. The input directory must contain at least
                                          one set of appropriate files for the pipeline to initiate.

              --merge [path/to/merge/dir]     Specify the path to a directory containing a subset of input reads for higher
                                          precision analysis. The pipeline will analyse the --input reads using 'erne-bs5'
                                          and the --merge reads using 'segemehl'. The file names in each directory must
                                          be identical, and all files must follow the same nomenclature as described for the
                                          --input parameter. Any number of read sets may be included, and they will be
                                          analysed in parallel. The individual reads specified in the --merge files must
                                          also be absent from the --input files. [default: off]

              --output [path/to/output/dir]   A path to a location to write the output results directory, which can be relative
                                          or absolute. This directory will contain sub-directories for each set of reads analysed
                                          during the pipeline. [default: wgbs]


         Options: REFERENCE GENOME [REQUIRED]
              --reference [path/to/genome.fa] Path to the input reference genome file in fasta format. EpiDiverse species are
                                          already available on the server in the '/scr/epi/genomes' directory. For species
                                          not within the scope of EpiDiverse, please refer to the 'Custom Reference Genome'
                                          section in the pipeline usage documentation.

              --thlaspi                       Shortcut to reference genome for 'Thlaspi arvense' on EpiDiverse infrastructure.
              --fragaria                      Shortcut to reference genome for 'Fragaria vesca' on EpiDiverse infrastructure.
              --populus                       Shortcut to reference genome for 'Populus nigra' on EpiDiverse infrastructure.


         Options: GENERAL

              --SE                            Interpret input reads as SINGLE END instead of PAIRED END. [default: off]

              --trim                          Enable read trimming step using 'cutadapt'. Additional parameters are specified
                                          below. [default: off]

              --fastqc                        Generate fastqc report for sequencing reads. Report will be generated for
                                          trimmed reads if the trimming process has been enabled. [default: off]

              --INDEX                         Specify if you would like the pipeline to generate the reference genome index
                                          automatically based on the options provided to the pipeline run. [default: off]

              --segemehl                      Enable bisulfite read alignment using only 'segemehl'. This has higher precision
                                          but is more memory and time intensive than 'erne-bs5'. [default: off]

              --unique                        Exclude multi-mapping reads from all post-alignment processing steps and the
                                          methylation extraction. [default: off]
                        
              --chrom [STR]                   Specify the scaffold/contig from which to estimate bisulfite non-conversion rate
                                          following alignment. Can be specified alongside eg. lambda, as the calculation will be
                                          made from the given contig without splitting the alignment file. [default: off]

              --split [path/to/fasta.fa]      If you have used a custom spike-in in place of 'lambda phage', then you must
                                          specify the path to the standalone fasta before splitting the alignments. NOTE: there
                                          must also be a fasta.fa.fai index file in the same directory. [default: lambda]
                                          
              --noLambda                      Disable the default processing/splitting of output files based on the assumption
                                          that 'lambda phage' DNA has been included during library preparation. This parameter
                                          dictates whether or not the reads mapping to the scaffold/contig specified in --split
                                          should be removed from the final alignment bam files. [default: off]
                  
              --noDedup                       Skip de-duplication step for downstream filtering of PCR duplicates. [default: off]

              --keepIndex                     Specify if you would like to keep generated index files. Can be used interchangeably
                                          with --INDEX to perform the same function while keeping index files. [default: off]

              --keepReads                     Specify if you would like to keep processed reads eg. after trimming [default: off]

              --keepBams                      Specify if you would like to keep intermediate bam files eg. raw alignments [default: off]


         Options: READ TRIMMING
              --forward [STR]                 The forward adapter sequence. [default: AGATCGGAAGAGCACACGTCTGAAC]
              --reverse [STR]                 The reverse adapter sequence. [default: AGATCGGAAGAGCGTCGTGTAGGGA]
              --clip5 [INT]                   Number of bases for hard clipping at the 5'-end. [default: 0]
              --clip3 [INT]                   Number of bases for hard clipping at the 3'-end. [default: 0]
              --minQual [INT]                 Minimum quality threshold for trimming. [default: 20]
              --minLeng [INT]                 Minimum length threshold of trimmed sequences. [default: 36]
              --minOver [INT]                 Minimum overlap length for adapter sequences in reads. [default: 3]


         Options: READ ALIGNMENT
              --minIns [INT]                 The minimum insert size between PE reads. [default: 0]
              --maxIns [INT]                 The maximum insert size between PE reads. [default: 500]

              --maxErrors [INT]                  Specify the maximum number of base errors allowed during read alignment with
                                          'erne-bs5'. This will disable the default auto-error function. [default: off]

              --minAccuracy [INT]                Specify the minimum percentage threshold for base matches during read alignment
                                          with 'segemehl'. [default: 90]

              --XF [FLOAT]                    Specify a proportion of opposite-strand bisulfite mismatches to allow before
                                          filtering. This excludes reads in regions that appear to have undergone genome
                                          rearrangement to the opposite strand, as alignment in this case is unreliable.
                                          [default: 0.03]


         Options: METHYLATION CONTEXT
              --noCpG                         Disable methylation calling in CpG context. [default: off]
              --noCHH                         Disable methylation calling in CHH context. [default: off]
              --noCHG                         Disable methylation calling in CHG context. [default: off]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode


         Example: 
              nextflow run epidiverse/wgbs \
              --input path/to/reads/dir \
              --reference path/to/reference.fa \
              --output wgbs \
              --unique

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         ===============================================
         E P I D I V E R S E - W G B S   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// DECLARE INITIAL PATH VARIABLES
// attempt to call check_ref_errors function from lib/functions.nf if workflow.profile is epidiverse
if ( workflow.profile.tokenize(",").contains("epi") || workflow.profile.tokenize(",").contains("diverse") ){

    // VALIDATE INITIAL PARAMETERS
    ParameterChecks.checkParamsEpi(params)
    include { check_ref_errors } from './lib/functions.nf' params(reference: params.reference, thlaspi: params.thlaspi, populus: params.populus, fragaria: params.fragaria, nolambda: params.noLambda)
    (fasta, fai, ebm_path, ctidx_path, gaidx_path) = check_ref_errors(params.reference, params.thlaspi, params.fragaria, params.populus, params.noLambda)
}

else {

    // VALIDATE INITIAL PARAMETERS
    ParameterChecks.checkParamsLocal(params)
    fasta = file("${params.reference}", checkIfExists: true, glob: false)
    fai = file("${params.reference}.fai", checkIfExists: true, glob: false)
    ebm_path = params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "${params.reference}/index/*.ebm" : "${params.reference}/lambda/*.ebm"
    ctidx_path = params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "${params.reference}/index/*.ctidx" : "${params.reference}/lambda/*.ctidx"
    gaidx_path = params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "${params.reference}/index/*.gaidx" : "${params.reference}/lambda/*.gaidx"
}

// VALIDATE ALL PARAMETERS
ParameterChecks.checkParams(params)

// establish path to reads in input and merge dirs
reads_path = params.SE ? "${params.input}/*.${params.extension}" : "${params.input}/*{1,2}.${params.extension}"
merge_path = params.SE ? "${params.merge}/*.${params.extension}" : "${params.merge}/*{1,2}.${params.extension}"
bam_path = "${params.input}/*.bam"

// determine contexts
if ((params.noCpG == true) && (params.noCHH == true) && (params.noCHG == true)) {error "ERROR: please specify methylation context for analysis"}
else if ((params.noCpG == true) && (params.noCHH == true) && (params.noCHG == false)) {context = "--noCpG --CHG "}
else if ((params.noCpG == true) && (params.noCHH == false) && (params.noCHG == true)) {context = "--noCpG --CHH "}
else if ((params.noCpG == true) && (params.noCHH == false) && (params.noCHG == false)) {context = "--noCpG --CHH --CHG "}
else if ((params.noCpG == false) && (params.noCHH == true) && (params.noCHG == true)) {context = " "}
else if ((params.noCpG == false) && (params.noCHH == true) && (params.noCHG == false)) {context = "--CHG "}
else if ((params.noCpG == false) && (params.noCHH == false) && (params.noCHG == true)) {context = "--CHH "}
else {context = "--CHH --CHG "}


// STAGE REFERENCE FILES
lamfa = file("${params.split}", checkIfExists: true, glob: false)
lai = file("${params.split}.fai", checkIfExists: true, glob: false)

ebm = params.INDEX || params.CALL ? Channel.empty() : file("${ebm_path}", checkIfExists: true)
ctidx = params.INDEX || params.CALL ? Channel.empty() : file("${ctidx_path}", checkIfExists: true)
gaidx = params.INDEX || params.CALL ? Channel.empty() : file("${gaidx_path}", checkIfExists: true)

// determine scaffold name of spike-in
chrom = lamfa.withReader{ it.readLine() }.tokenize(' ').get(0).substring(1)

// PRINT LOGGING INFO
if (params.CALL && !params.WGBS){

    // PRINT SECONDARY LOGGING INFO
    log.info ""
    log.info "         ================================================="
    log.info "          E P I D I V E R S E - C A L L   P I P E L I N E"
    if (params.debug){
    log.info "         (debug mode enabled)"
    log.info "         =================================================" }
    else {
    log.info "         =================================================" }
    log.info "         ~ version ${workflow.manifest.version}"
    log.info ""
    log.info "         input dir      : ${params.input}"
    log.info "         output dir     : ${params.output}"
    log.info "         PCR dups       : ${params.noDedup ? "ignored" : "marked" }"
    log.info "         context(s)     : ${params.noCpG ? "" : "CpG " }${params.noCHH ? "" : "CHH " }${params.noCHG ? "" : "CHG" }"
    log.info "         chrom          : ${params.chrom ? "${params.chrom} " : "-" }"
    log.info ""

} else {

    // PRINT PRIMARY LOGGING INFO
    log.info ""
    log.info "         ================================================="
    log.info "          E P I D I V E R S E - W G B S   P I P E L I N E"
    if (params.debug){
        log.info "         (debug mode enabled)"
        log.info "         ================================================="
    } else {
        log.info "         ================================================="
    }
    log.info "         ~ version ${workflow.manifest.version}"
    log.info ""
    log.info "         reference      : ${fasta.baseName}"
    log.info "         input dir      : ${workflow.profile.contains("test") ? "test profile" : "${params.input}"}"
    if (params.merge){
        log.info "         merge dir      : ${workflow.profile.contains("test") ? "test profile" : "${params.merge}"}"
    }
    log.info "         output dir     : ${params.output}"
    log.info "         extension      : *.${params.extension}"
    log.info "         read type      : ${params.SE ? "single-end" : "paired-end (min: $params.minIns max: $params.maxIns)" }"
    log.info "         read trimming  : ${params.trim ? 'enabled' : 'disabled' }"
    log.info "         fastqc report  : ${params.fastqc ? 'enabled' : 'disabled' }"
    log.info "         mapping tool   : ${params.merge ? "combined" : params.segemehl ? "segemehl" : "erne-bs5"}"
    log.info "         multimappings  : ${params.unique ? "filtered" : "included" }"
    log.info "         PCR dups       : ${params.noDedup ? "included" : "filtered" }"
    log.info "         XF filter      : ${params.XF}"
    log.info "         conv.rate est. : ${params.chrom ? "${params.chrom} " : "" }${params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "" : "${chrom}" }"
    log.info "         context(s)     : ${params.noCpG ? "" : "CpG " }${params.noCHH ? "" : "CHH " }${params.noCHG ? "" : "CHG" }"
    log.info ""

    // PRINT ADDITIONAL LOGGING INFO
    if (params.trim){
    log.info "         ================================================="
    log.info "         trimming options"
    log.info "         ================================================="
    log.info "         forward adptr  : ${params.forward}"
    log.info "         reverse adptr  : ${params.reverse}"
    log.info "         hard clip 5\'   : ${!params.clip5 ? "disabled" : "$params.clip5 bases" }"
    log.info "         hard clip 3\'   : ${!params.clip3 ? 'disabled' : "$params.clip3 bases" }"
    log.info "         min. quality   : ${params.minQual}"
    log.info "         min. length    : ${params.minLeng}"
    log.info "         min. overlap   : ${params.minOver}"
    log.info ""
    }
}

log.info "         ================================================"
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



/////////////////////
// COMMON CHANNELS //
/////////////////////

// CALL workflow takes priority over -profile test
if ( params.CALL && !params.WGBS ){

    // STAGE BAM CHANNEL
    bam = Channel
        .fromFilePairs(bam_path, size: 1)
        .ifEmpty{ exit 1, "ERROR: cannot find valid BAM files in dir: ${params.input}\n \
        Did you mean to specify the --CALL option?"}
        .map {tuple(it[0], "subset", it[1].flatten())}

} else {

    // attempt to call check_test_data function from lib/functions.nf if workflow.profile contains test
    if ( workflow.profile.tokenize(",").contains("test") ){

        include { check_test_data } from './lib/functions.nf' params(readPaths: params.readPaths, mergePaths: params.mergePaths, singleEnd: params.SE, merge: params.merge)
        (reads, merged) = check_test_data(params.readPaths, params.mergePaths, params.SE, params.merge)
        
    } else { 

        // STAGE READS INPUT CHANNEL
        reads = Channel
            .fromFilePairs(reads_path, size: params.SE ? 1 : 2)
            .ifEmpty{ exit 1, "ERROR: cannot find valid read files in dir: ${params.input}\n \
            The pipeline will expect PE reads in compressed *{1,2}.${params.extension} format\n \
            unless you have specified the --SE parameter or a different extension using --extension"}
            .map{ tuple(it[0], "input", it[1]) }
            .take(params.take.toInteger())

        //STAGE READS MERGE CHANNEL
        merged = !params.merge ? Channel.empty() :
        Channel
            .fromFilePairs(merge_path, size: params.SE ? 1 : 2)
            .ifEmpty{ exit 1, "ERROR: cannot find valid read files in dir: ${params.merge}\n \
            The pipeline will expect PE reads in compressed *{1,2}.${params.extension} format\n \
            unless you have specified the --SE parameter or a different extension using --extension"}
            .map{ tuple(it[0], "merge", it[1]) }
            .take(params.take.toInteger())
    }
}


////////////////////
// BEGIN PIPELINE //
////////////////////


// INCLUDES
include {erne_bs5_indexing;segemehl_indexing} from './lib/index.nf' params(params)
include {read_trimming;read_merging;fastqc;erne_bs5;segemehl;erne_bs5_processing;segemehl_processing;bam_merging;bam_subsetting;bam_statistics;bam_filtering;} from './lib/wgbs.nf' params(params)
//include {bam_grouping;bam_sampling;bam_processing;Picard_MarkDuplicates;MethylDackel;conversion_rate_estimation} from './lib/call.nf' params(params)
include {bam_processing;Picard_MarkDuplicates;MethylDackel;conversion_rate_estimation} from './lib/call.nf' params(params)


// WORKFLOWS

// INDEX workflow - secondary pipeline
workflow 'INDEX' {

    take:
        fasta
        fai
        lamfa
        lai

    main:
        erne_bs5_indexing(fasta,fai,lamfa,lai)
        segemehl_indexing(fasta,fai,lamfa,lai)

    emit:
        ebm = erne_bs5_indexing.out
        ctidx = segemehl_indexing.out[0]
        gaidx = segemehl_indexing.out[1]
}

// WGBS workflow - primary pipeline
workflow 'WGBS' {

    take:
        reads
        merged
        ebm
        ctidx
        gaidx
        fasta
        fai
        lamfa
        lai
        chrom
 
    main:
        // read trimming and merging
        read_trimming(reads.mix(merged))
        params.trim ? read_merging(read_trimming.out[0].groupTuple().map{ tuple(it[0], it[1], *it[2]) }) :\
        read_merging(reads.mix(merged).groupTuple().map{ tuple(it[0], it[1], *it[2]) })

        // fastqc process
        params.merge ? fastqc(read_merging.out[0]) :\
        params.trim ? fastqc(read_trimming.out[0]) : fastqc(reads)

        // erne_bs5 process
        params.trim ? erne_bs5(read_trimming.out[0], ebm) : erne_bs5(reads, ebm)

        // segemehl process
        params.trim ? segemehl(read_trimming.out[0], ctidx, gaidx, fasta, lamfa) : params.merge ?\
        segemehl(merged, ctidx, gaidx, fasta, lamfa) : segemehl(reads, ctidx, gaidx, fasta, lamfa)

        // alignment post-processing
        erne_bs5_processing(erne_bs5.out[0],fasta,lamfa)
        segemehl_processing(segemehl.out[0])

        // alignment merging and subsetting
        bam_merging(erne_bs5_processing.out[0].combine(segemehl_processing.out[0], by: 0))
        params.merge ? bam_subsetting(bam_merging.out[0],fai,lai,chrom) : bam_subsetting(erne_bs5_processing.out[0].mix(segemehl_processing.out[0]),fai,lai,chrom)
        !params.noLambda || params.split != "${baseDir}/data/lambda.fa" ? bam_filtering(bam_subsetting.out[0].mix(bam_subsetting.out[1])) : params.merge ?\
        bam_filtering(bam_merging.out[0]) : bam_filtering(erne_bs5_processing.out[0].mix(segemehl_processing.out[0]))

        // alignment statistics
        !params.noLambda || params.split != "${baseDir}/data/lambda.fa" ? bam_statistics(bam_subsetting.out[1]) : params.merge ?\
        bam_statistics(bam_merging.out[0]) : bam_statistics(erne_bs5_processing.out[0].mix(segemehl_processing.out[0]))

    emit:
        erne_bs5_processing_publish = erne_bs5_processing.out[0]
        segemehl_processing_publish = segemehl_processing.out[0]
        bam_merging_publish = bam_merging.out[0]
        bam_subsetting_publish_lambda = bam_subsetting.out[0]
        bam_subsetting_publish_subset = bam_subsetting.out[1]
        bam_filtering_out = bam_filtering.out[0]

}

// CALL workflow - secondary pipeline for Methylation Quantification
workflow "CALL" {

    take:
        bam
        fasta
        lamfa
        context
        chrom

    main:
        /*
        bam1 = bam.filter{ it[1] != "lambda" }
        bam2 = bam.map{ tuple(*it, it[0]) }
        // split by read groups
        bam_grouping(bam1)
        //bam_grouping_out = bam_grouping.out.map{ tuple(it[0], tuple(it[1].flatten())) }
        txt1 = bam_grouping.out.filter{ it[1] instanceof Path }.map{ tuple(it[0], it[0]) } // only one RG present
        txt2 = bam_grouping.out.filter{ it[1] instanceof Collection && it[1].size() > 1  }.transpose().map{ tuple(it[0], it[1].baseName.tokenize(".").init().join(""), it[1]) }
        rgs1 = bam.filter{ it[1] != "lambda" }.combine(txt1, by: 0) // eg. [replicate, bamtype, *.bam, filename]
        rgs2 = bam.filter{ it[1] != "lambda" }.combine(txt2, by: 0) // eg. [replicate, bamtype, *.bam, filename, *.txt]
        bam_sampling(rgs2) // eg. [replicate, bamtype, sample.bam, filename]
       
        // deduplication and methylation calling
        params.splitRG ? bam_processing(bam2.filter{ it[1] == "lambda" }.mix(rgs1, bam_sampling.out)) : bam_processing(bam2)
        Picard_MarkDuplicates(bam_processing.out)
        params.noDedup ? MethylDackel(bam_processing.out,fasta,lamfa,context) : MethylDackel(Picard_MarkDuplicates.out[0],fasta,lamfa,context)

        // conversion rate estimation and duplication statistics
        conversion_rate_estimation(MethylDackel.out[0],chrom)
        */

        // deduplication and methylation calling
        bam_processing(bam)
        Picard_MarkDuplicates(bam_processing.out)
        params.noDedup ? MethylDackel(bam_processing.out,fasta,lamfa,context) : MethylDackel(Picard_MarkDuplicates.out[0],fasta,lamfa,context)

        // conversion rate estimation and duplication statistics
        conversion_rate_estimation(MethylDackel.out[0],chrom)



    emit:
        conversion_rate_publish = conversion_rate_estimation.out
    
}


// MAIN workflow
workflow {

    main:
        // skip INDEX and WGBS workflow
        if (params.CALL && !params.WGBS) {

            CALL(bam,fasta,lamfa,context,chrom)
            CALL.out.conversion_rate_publish.collectFile().subscribe{ it.copyTo("${params.output}/bam/${it.baseName}/stats/BisNonConvRate.txt") }

            //INDEX(Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty())
            //WGBS(Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty())

        } else {

            // genome index required
            if (params.INDEX) {

                INDEX(fasta,fai,lamfa,lai)
                WGBS(reads,merged,INDEX.out.ebm,INDEX.out.ctidx,INDEX.out.gaidx,fasta,fai,lamfa,lai,chrom)

            // genome index available
            } else {

                //INDEX(Channel.empty(),Channel.empty(),Channel.empty(),Channel.empty())
                WGBS(reads,merged,ebm,ctidx,gaidx,fasta,fai,lamfa,lai,chrom)

            }

            // filter unique alignments
            if (params.unique) {
            
                bam = WGBS.out.bam_filtering_out

            // no filtering
            } else {

                bam = !params.noLambda || params.split != "${baseDir}/data/lambda.fa" ? WGBS.out.bam_subsetting_publish_lambda.mix(WGBS.out.bam_subsetting_publish_subset) : params.merge ?\
                WGBS.out.bam_merging_publish : WGBS.out.erne_bs5_processing_publish.mix(WGBS.out.segemehl_processing_publish)

            }

            CALL(bam,fasta,lamfa,context,chrom)
            CALL.out.conversion_rate_publish.collectFile().subscribe{ it.copyTo("${params.output}/bam/${it.baseName}/stats/BisNonConvRate.txt") }

        }

        /*
        // CALL workflow
        if (params.CALL || !params.WGBS) {
            CALL(bam,fasta,lamfa,context,chrom)
            CALL.out.conversion_rate_publish.collectFile().subscribe{ it.copyTo("${params.output}/bam/${it.baseName}/stats/BisNonConvRate.txt") }
        }
        */

}


//////////////////
// END PIPELINE //
//////////////////


// WORKFLOW TRACING
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${workflow.success && !params.debug ? "(cleared)" : ""}"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    if (workflow.success && !params.debug) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}