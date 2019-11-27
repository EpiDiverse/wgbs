#!/usr/bin/env nextflow

// initial vars
cutadapt_clip5 = params.clip5.toInteger() > 0 ? " -u ${params.clip5}" : ""
cutadapt_clip3 = params.clip3.toInteger() > 0 ? " -u -${params.clip3}" : ""
erne_errors = params.maxErrors.toInteger() < 0 ? "--errors ${params.maxErrors} " : " "



// STAGE INPUT READS INTO DIRECTORIES
// requires reads from params.input
process "stage_input_directories" {

    label 'tiny'
    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, path(reads)

    output:
    tuple replicate, val("input"), path("$replicate")

    script:
    """
    mkdir ${replicate} ${replicate}/fastq
    cp -a *.${params.extension} ${replicate}/fastq
    """

}



// STAGE MERGE READS INTO DIRECTORIES
// requires reads from params.merge
process "stage_merge_directories" {

    label 'tiny'
    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, path(reads)

    output:
    tuple replicate, val("merge"), path("$replicate")

    when:
    params.merge

    script:
    """
    mkdir ${replicate} ${replicate}/fastq
    cp -a *.${params.extension} ${replicate}/fastq
    """

}



// TRIM READS WITH CUTADAPT (OPTIONAL)
// requires build_input_directories.mix(build_merge_directories)
process "read_trimming" {

    label 'low'
    label 'finish'
    tag "$replicate"
    
    //publishDir "$results_path/$replicate/fastqc", saveAs: { params.merge == params.input ? "${it}" : null }, mode: 'copy'

    input:
    tuple replicate, readtype, path("inputs")
    // eg. [replicate, [read1.fastq.gz, read2.fastq.gz]]

    output:
    tuple replicate, readtype, path("$replicate")
    // eg. [replicate, /path/to/replicate]
    path "$replicate/fastq/*.${params.extension}"
    path "$replicate/fastq/logs/*.log"

    when:
    params.trim

    script:
    if( params.SE )
        """
        fq=inputs/fastq/*.${params.extension}
        mkdir ${replicate} ${replicate}/fastq ${replicate}/fastq/logs

        cutadapt -j ${task.cpus} -a ${params.forward}${cutadapt_clip5}${cutadapt_clip3} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o ${replicate}/fastq/\$(basename \$fq) \$fq \\
        > ${replicate}/fastq/logs/${readtype}.log 2>&1
        """
    else
        """
        fq1=inputs/fastq/*1.${params.extension}
        fq2=inputs/fastq/*2.${params.extension}
        mkdir ${replicate} ${replicate}/fastq ${replicate}/fastq/logs

        cutadapt -j ${task.cpus} -a ${params.forward} -A ${params.reverse}${cutadapt_clip5}${cutadapt_clip3} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o ${replicate}/fastq/\$(basename \$fq1) -p ${replicate}/fastq/\$(basename \$fq2) \$fq1 \$fq2 \\
        > ${replicate}/fastq/logs/${readtype}.log 2>&1
        """

}



// COMBINE INPUT AND MERGE READS FOR FASTQC AND/OR PUBLISHING AFTER TRIMMING
// requires groupTuple() by replicate 
process "read_merging" {

    label 'tiny'
    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, readtype, path("inputs*")
    // eg. [replicate, ["input","merge"], [/path/to/input/replicate, /path/to/merge/replicate]]

    output:
    tuple replicate, readtype, path("$replicate")
    // eg. [replicate, ["input","merge"], /path/to/replicate]
    path "$replicate/fastq/*.${params.extension}"

    when:
    (params.merge && (params.fastqc || (params.trim && params.keepReads)))

    script:
    if( params.SE )
        """
        mkdir ${replicate} ${replicate}/fastq
        cat inputs*/fastq/*.${params.extension} > ${replicate}/fastq/${replicate}.${params.extension}
        """
    else
        """
        mkdir ${replicate} ${replicate}/fastq
        echo -e "1\\n2" | xargs -I{} sh -c 'cat inputs*/fastq/*\$1.${params.extension} > ${replicate}/fastq/${replicate}_\$1.${params.extension}' -- {}
        """
}



// VIEW READS QC WITH FASTQC (OPTIONAL)
process "fastqc" {

    label 'tiny'
    label 'ignore'
    tag "$replicate"

    input:
    tuple replicate, readtype, path("inputs")
    // eg. [replicate, "input", /path/to/replicate]
    // eg. [replicate, ["input","merge"], /path/to/replicate]

    output:
    path "$replicate/fastq/*"
    path "$replicate/fastq/logs/fastqc.log"

    when:
    params.fastqc

    script:
    """
    mkdir ${replicate} ${replicate}/fastq ${replicate}/fastq/logs
    fastqc inputs/fastq/*.${params.extension} -threads ${task.cpus} -outdir=${replicate}/fastq > ${replicate}/fastq/logs/fastqc.log
    """

}



// FAST READ ALIGNMENT USING "erne-bs5" 
process "erne_bs5" {

    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, readtype, path("inputs")
    // eg. [replicate, ["input"], /path/to/inputs]
    path ebm

    output:
    tuple replicate, path("$replicate/bam/raw.erne-bs5.bam")
    // eg. [replicate, /path/to/replicate/bam/raw.erne-bs5.bam]
    //path "$replicate/bam/raw.erne-bs5.bam"
    path "$replicate/bam/logs/raw.erne-bs5.log"

    when:
    ( !params.segemehl || params.merge ) && readtype == "input"

    script:
    if( params.SE )
        """
        fq=inputs/fastq/*.${params.extension}
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs

        erne-bs5 --reference ${ebm} --query1 \$fq --fragment-size-min ${params.minIns} --fragment-size-max ${params.maxIns} \\
        ${erne_errors}--threads ${task.cpus - 2} --output unsorted.erne-bs5.bam --print-all \\
        > ${replicate}/bam/logs/raw.erne-bs5.log 2>&1 || exit \$?
        samtools sort -o ${replicate}/bam/raw.erne-bs5.bam unsorted.erne-bs5.bam
        """
    else
        """
        fq1=inputs/fastq/*1.${params.extension}
        fq2=inputs/fastq/*2.${params.extension}
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs

        erne-bs5 --reference ${ebm} --query1 \$fq1 --query2 \$fq2 --fragment-size-min ${params.minIns} --fragment-size-max ${params.maxIns} \\
        ${erne_errors}--threads ${task.cpus - 2} --output unsorted.erne-bs5.bam --print-all \\
        > ${replicate}/bam/logs/raw.erne-bs5.log 2>&1 || exit \$?
        samtools sort -o ${replicate}/bam/raw.erne-bs5.bam unsorted.erne-bs5.bam
        """

}



// SENSITIVE READ ALIGNMENT USING "segemehl"
process "segemehl" {

    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, readtype, path("inputs")
    // eg. [replicate, ["input"], /path/to/inputs]
    path ctidx
    path gaidx
    path fasta
    path lamfa

    output:
    tuple replicate, path("$replicate/bam/raw.segemehl.bam")
    // eg. [replicate, /path/to/replicate]
    //path "$replicate/bam/raw.segemehl.bam"
    path "$replicate/bam/logs/raw.segemehl.log" 

    when:
    params.segemehl || ( params.merge && readtype == "merge" )

    script:
    if( params.SE )
        """
        fq=inputs/fastq/*.${params.extension}
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs

        segemehl.x -i ${ctidx} -j ${gaidx} \\
        -d ${params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "${fasta}" : "${fasta} ${lamfa}"} \\
        -q \$fq -o raw.segemehl.sam -I ${params.maxIns} -A ${params.minAccuracy} -s -t ${task.cpus} -F 1 -H 1 -D 1 \\
        > ${replicate}/bam/logs/raw.segemehl.log 2>&1 || exit \$?
        samtools view -Sb raw.segemehl.sam > ${replicate}/bam/raw.segemehl.bam
        """
    else
        """
        fq1=inputs/fastq/*1.${params.extension}
        fq2=inputs/fastq/*2.${params.extension}
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs
        
        segemehl.x -i ${ctidx} -j ${gaidx} \\
        -d ${params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "${fasta}" : "${fasta} ${lamfa}"} \\
        -q \$fq1 -p \$fq2 -o raw.segemehl.sam -I ${params.maxIns} -A ${params.minAccuracy} -s -t ${task.cpus} -F 1 -H 1 -D 1 \\
        > ${replicate}/bam/logs/raw.segemehl.log 2>&1 || exit \$?
        samtools view -Sb raw.segemehl.sam > ${replicate}/bam/raw.segemehl.bam
        """

}



// POST-PROCESS BAMFILES READY FOR DOWNSTREAM ANALYSIS (erne-bs5)
process "erne_bs5_processing" {

    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, path(erne)
    // eg. [replicate, /path/to/input]
    path fasta
    path lamfa
    
    output:
    tuple replicate, val("erne-bs5"), path("$replicate/bam/proc.erne-bs5.bam")
    // eg. [replicate, /path/to/replicate]
    path "${replicate}/${replicate}.bam"
    
    script:
    if (params.SE && params.noLambda && params.split == "${baseDir}/data/lambda.fa")
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s bam/proc.erne-bs5.bam ${replicate}/${replicate}.bam

        samtools index ${erne}
        correct_sam_cigar.py ${erne} corrected.erne-bs5.bam || exit \$?
    
        samtools sort -T deleteme -o sorted.erne-bs5.bam corrected.erne-bs5.bam
        samtools index sorted.erne-bs5.bam
        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . ${fasta} sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
    else if (!params.SE && params.noLambda && params.split == "${baseDir}/data/lambda.fa")
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s bam/proc.erne-bs5.bam ${replicate}/${replicate}.bam

        samtools index ${erne}
        correct_sam_format.py -i ${params.maxIns} -t ${task.cpus} -T . ${erne} corrected.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -no unsorted.erne-bs5.bam corrected.erne-bs5.bam
        correct_sam_tlens -i unsorted.erne-bs5.bam > tlens.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -o sorted.erne-bs5.bam tlens.erne-bs5.bam
        samtools index sorted.erne-bs5.bam
        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . ${fasta} sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
    else if (params.SE)
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s bam/proc.erne-bs5.bam ${replicate}/${replicate}.bam

        if [ -z "\$(tail -c 1 ${fasta})" ]
        then
        cat ${fasta} ${lamfa} > fasta.tmp
        else
        cat ${fasta} > fasta.tmp; echo >> fasta.tmp; cat ${lamfa} >> fasta.tmp
        fi

        samtools index ${erne}
        correct_sam_cigar.py ${erne} corrected.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -o sorted.erne-bs5.bam corrected.erne-bs5.bam
        samtools index sorted.erne-bs5.bam
        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . fasta.tmp sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
    else
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s bam/proc.erne-bs5.bam ${replicate}/${replicate}.bam

        if [ -z "\$(tail -c 1 ${fasta})" ]
        then
        cat ${fasta} ${lamfa} > fasta.tmp  
        else
        cat ${fasta} > fasta.tmp; echo >> fasta.tmp; cat ${lamfa} >> fasta.tmp      
        fi

        samtools index ${erne}
        correct_sam_format.py -i ${params.maxIns} -t ${task.cpus} -T . ${erne} corrected.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -no unsorted.erne-bs5.bam corrected.erne-bs5.bam
        correct_sam_tlens -i unsorted.erne-bs5.bam > tlens.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -o sorted.erne-bs5.bam tlens.erne-bs5.bam
        samtools index sorted.erne-bs5.bam

        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . fasta.tmp sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
}



// POST-PROCESS BAMFILES READY FOR DOWNSTREAM ANALYSIS (segemehl)
process "segemehl_processing" {

    label 'low'
    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, path(sege)
    // eg. [replicate, /path/to/input]

    output:
    tuple replicate, val("segemehl"), path("$replicate/bam/proc.segemehl.bam")
    // eg. [replicate, segemehl, /path/to/proc.segemehl.bam]
    path "${replicate}/${replicate}.bam"

    script:
    if(params.SE)
        """
        mkdir ${replicate} ${replicate}/bam
        samtools sort -o sorted.segemehl.bam ${sege}
        samtools index sorted.segemehl.bam

        filter_sam_xf_tag.py -c ${params.XF} sorted.segemehl.bam ${replicate}/bam/proc.segemehl.bam || exit \$?
        ln -s bam/proc.segemehl.bam ${replicate}/${replicate}.bam
        """
    else
        """
        mkdir ${replicate} ${replicate}/bam
        correct_sam_tlens -i ${sege} > tlens.segemehl.bam || exit \$?
        samtools sort -o sorted.segemehl.bam tlens.segemehl.bam
        samtools index sorted.segemehl.bam
        
        filter_sam_xf_tag.py -c ${params.XF} sorted.segemehl.bam ${replicate}/bam/proc.segemehl.bam || exit \$?
        ln -s bam/proc.segemehl.bam ${replicate}/${replicate}.bam
        """
}



// MERGE BAMFILES READY FOR DOWNSTREAM ANALYSIS (OPTIONAL)
process "bam_merging" {

    label 'low'
    label 'finish'
    tag "$replicate"
    
    input:
    tuple replicate, erne, path(erne_bs5), sege, path(segemehl)
    // eg. [replicate, proc.erne-bs5.bam, proc.segemehl.bam]

    output:
    tuple replicate, val("merged"), path("$replicate/bam/merged.bam")
    // eg. [replicate, merged, /path/to/replicate/replicate.bam]
    path "${replicate}/${replicate}.bam"
   
    when:
    params.merge
    
    script:
    """
    mkdir ${replicate} ${replicate}/bam
    ln -s bam/merged.bam ${replicate}/${replicate}.bam
    samtools merge ${replicate}/bam/merged.bam ${erne_bs5} ${segemehl}
    """
}



// SPLIT BAMFILES FOR LAMBDA AND EVERTHING ELSE
process "bam_subsetting" {

    label 'low'
    label 'finish'
    tag "$replicate"

    input:
    tuple replicate, bamtype, path(bamfile)
    // eg. [replicate, merged, bamfile.bam]
    path fai
    path lai
    val chrom

    output:
    tuple replicate, val("lambda"), path("${replicate}/bam/${chrom}.bam")
    tuple replicate, val("subset"), path("${replicate}/bam/subset.bam")
    // eg. [replicate, subset, /path/to/replicate/bam/subset.bam]
    path "${replicate}/${replicate}.bam"
    
    when:
    !params.noLambda || params.split != "${baseDir}/data/lambda.fa"
   
    script:
    """
    mkdir ${replicate} ${replicate}/bam
    ln -s bam/subset.bam ${replicate}/${replicate}.bam

    # reset alignments whose mates have erroneously aligned to lambda
    samtools view -h ${bamfile} |
    awk 'BEGIN {OFS="\\t"} {if((\$3=="${chrom}" && \$7=="=") || \$7=="${chrom}") \
    {if(\$2==97){\$2=73} else { if(\$2==81){\$2=89} \
    else { if(\$2==161){\$2=137} else{ if(\$2==145){\$2=153} }}}; \$7="*"; \$8=0}; {print \$0}}' |
    samtools sort -T deleteme -o sort.bam - || exit \$?
    samtools index sort.bam

    # remove lambda from the header and split bams
    samtools view -H sort.bam | awk '\$1~/^@[[:upper:]]{2}/ && \$2!="SN:${chrom}"' > header.sam || exit \$?
    awk '{printf("%s\\t0\\t%s\\n",\$1,\$2);}' ${fai} > fai.tmp
    samtools view -bL fai.tmp sort.bam | samtools reheader header.sam - > ${replicate}/bam/subset.bam || exit \$?
    samtools view sort.bam ${chrom} | samtools view -bt ${lai} - > ${replicate}/bam/${chrom}.bam
    """
}



// CALCULATE ALIGNMENT STATISTICS
process "bam_statistics" {

    label 'low'
    label 'ignore'
    tag "$replicate"
    
    input:
    tuple replicate, bamtype, path(bamfile)
    // eg. [replicate, bamtype, bamfile.bam]

    output:
    path "${replicate}/stats/${replicate}.bam.stats"
    path "${replicate}/stats/bam/*.png"    

    script:
    """
    mkdir ${replicate} ${replicate}/stats ${replicate}/stats/bam
    samtools sort -T deleteme -o sorted.bam ${bamfile}
    samtools stats sorted.bam > ${replicate}/stats/${replicate}.bam.stats
    plot-bamstats -p ${replicate}/stats/bam/ ${replicate}/stats/${replicate}.bam.stats
    """    
}



// PRE-PROCESS BAMFILES READY FOR DOWNSTREAM METHYLATION EXTRACTION
process "bam_processing" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"
    
    input:
    tuple replicate, bamtype, path(bamfile)
    // eg. [replicate, lambda, /path/to/bamfile.bam]
    // eg. [replicate, subset, /path/to/bamfile.bam]

    output:
    tuple replicate, bamtype, path("$replicate/bam/unique.bam")
    tuple replicate, bamtype, path("$replicate/bam/*.bam")
    // eg. [replicate, lambda, /path/to/replicate/*.bam]
    // eg. [replicate, subset, /path/to/replicate/*.bam]
    
    script:
    if (params.unique)
        """
        mkdir ${replicate} ${replicate}/bam
        change_sam_qname -i ${bamfile} -o ${replicate}/bam/unique.bam --tags HI XB --read_name_tag XN
        """
    else
        """
        mkdir ${replicate} ${replicate}/bam
        filter_sam_uniqs.py ${bamfile} filtered.bam ${replicate}/bam/multimapped.bam
        change_sam_qname -i filtered.bam -o ${replicate}/bam/unique.bam --tags HI XB --read_name_tag XN
        """
}
