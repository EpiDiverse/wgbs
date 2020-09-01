#!/usr/bin/env nextflow

// initial vars
//cutadapt_clip5 = params.clip5.toInteger() > 0 ? " -u ${params.clip5}" : ""
//cutadapt_clip3 = params.clip3.toInteger() > 0 ? " -u -${params.clip3}" : ""
//erne_errors = params.maxErrors.toInteger() < 0 ? "--errors ${params.maxErrors} " : " "



// TRIM READS WITH CUTADAPT (OPTIONAL)
// requires build_input_directories.mix(build_merge_directories)
process "read_trimming" {

    label 'low'
    label 'finish'
    tag "$replicate"

    publishDir "${params.output}", pattern: "fastq/*.${params.extension}", mode: 'copy', enabled: params.keepReads && !params.merge ? true : false
    publishDir "${params.output}", pattern: "fastq/logs/*.log", mode: 'move'

    input:
    tuple val(replicate), val(readtype), path(reads)
    // eg. [replicate, "input", [read1.fastq.gz, read2.fastq.gz]]

    output:
    tuple val(replicate), val(readtype), path("fastq/*.${params.extension}")
    // eg. [replicate, "input", /path/to/replicate]
    path "fastq/*.${params.extension}"
    path "fastq/logs/*.log"

    when:
    params.trim

    script:
    if( params.SE )
        """
        mkdir fastq fastq/logs
        cutadapt -j ${task.cpus} -a ${params.forward}${params.clip5.toInteger() > 0 ? " -u ${params.clip5}" : ""}${params.clip3.toInteger() > 0 ? " -u -${params.clip3}" : ""} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o fastq/${params.merge ? "${readtype}." : ""}${replicate}.${params.extension} ${reads} \\
        > fastq/logs/cutadapt.${replicate}.${readtype}.log 2>&1
        """
    else
        """
        mkdir fastq fastq/logs
        cutadapt -j ${task.cpus} -a ${params.forward} -A ${params.reverse}${params.clip5.toInteger() > 0 ? " -u ${params.clip5}" : ""}${params.clip3.toInteger() > 0 ? " -u -${params.clip3}" : ""} \\
        -q ${params.minQual} -m ${params.minLeng} -O ${params.minOver} \\
        -o fastq/${params.merge ? "${readtype}." : ""}${reads[0]} \\
        -p fastq/${params.merge ? "${readtype}." : ""}${reads[1]} ${reads} \\
        > fastq/logs/cutadapt.${replicate}.${readtype}.log 2>&1
        """

}



// COMBINE INPUT AND MERGE READS FOR FASTQC AND/OR PUBLISHING AFTER TRIMMING
// requires groupTuple() by replicate 
process "read_merging" {

    label 'tiny'
    label 'finish'
    tag "$replicate"

    publishDir "${params.output}", mode: 'copy', enabled: params.keepReads && params.trim ? true : false

    input:
    tuple val(replicate), val(readtype), path("input"), path("merge")
    // eg. [replicate, ["input","merge"], [/path/to/input/replicate, /path/to/merge/replicate]]

    output:
    tuple val(replicate), val(readtype), path("fastq/*.${params.extension}")
    // eg. [replicate, ["input","merge"], /path/to/replicate]
    path "fastq/*.${params.extension}"

    when:
    (params.merge && (params.fastqc || (params.trim && params.keepReads)))

    script:
    if( params.SE )
        """
        mkdir fastq
        cat input* merge* > fastq/${replicate}.${params.extension}
        """
    else
        """
        mkdir fastq
        echo -e "1\\n2" | xargs -I{} sh -c 'cat *\$1 > fastq/${replicate}_\$1.${params.extension}' -- {}
        """
}



// VIEW READS QC WITH FASTQC (OPTIONAL)
process "fastqc" {

    label 'tiny'
    label 'ignore'
    tag "$replicate"

    publishDir "${params.output}", mode: 'move'

    input:
    tuple val(replicate), val(readtype), path(reads)
    // eg. [replicate, "input", /path/to/replicate]
    // eg. [replicate, ["input","merge"], /path/to/replicate]

    output:
    path "fastq/*.{html,zip}"
    path "fastq/logs/fastqc.${replicate}.log"

    when:
    params.fastqc

    script:
    """
    mkdir fastq fastq/logs
    fastqc ${reads} -threads ${task.cpus} -outdir=fastq > fastq/logs/fastqc.${replicate}.log
    """

}



// FAST READ ALIGNMENT USING "erne-bs5" 
process "erne_bs5" {

    label 'finish'
    tag "$replicate"

    publishDir "${params.output}/bam", pattern: "$replicate/bam/raw.erne-bs5.bam", mode: 'copy', enabled: params.keepBams ? true : false
    publishDir "${params.output}/bam", pattern: "$replicate/bam/logs/raw.erne-bs5.log", mode: 'move'

    input:
    tuple val(replicate), val(readtype), path(reads)
    // eg. [replicate, ["input"], /path/to/inputs]
    path ebm

    output:
    tuple val(replicate), path("$replicate/bam/raw.erne-bs5.bam")
    // eg. [replicate, /path/to/replicate/bam/raw.erne-bs5.bam]
    //path "$replicate/bam/raw.erne-bs5.bam"
    path "$replicate/bam/logs/raw.erne-bs5.log"

    when:
    ( !params.segemehl || params.merge ) && readtype == "input"

    script:
    if( params.SE )
        """
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs

        erne-bs5 --reference ${ebm} --query1 ${reads} --fragment-size-min ${params.minIns} --fragment-size-max ${params.maxIns} \\
        ${params.maxErrors.toInteger() < 0 ? "--errors ${params.maxErrors} " : " "}--threads ${task.cpus - 2} --output unsorted.erne-bs5.bam --print-all \\
        > ${replicate}/bam/logs/raw.erne-bs5.log 2>&1 || exit \$?

        samtools view -H unsorted.erne-bs5.bam > header.txt
        for ((i=1;i<=\$(grep -c SM:no_sample_specified <(samtools view -H unsorted.erne-bs5.bam));i++));
        do sed -i "0,/SM:no_sample_specified/{s/SM:no_sample_specified/SM:sample\${i}/}" header.txt; done || exit \$?

        samtools reheader header.txt unsorted.erne-bs5.bam |
        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o ${replicate}/bam/raw.erne-bs5.bam -
        """
    else
        """
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs

        erne-bs5 --reference ${ebm} --query1 ${reads[0]} --query2 ${reads[1]} --fragment-size-min ${params.minIns} --fragment-size-max ${params.maxIns} \\
        ${params.maxErrors.toInteger() < 0 ? "--errors ${params.maxErrors} " : " "}--threads ${task.cpus - 2} --output unsorted.erne-bs5.bam --print-all \\
        > ${replicate}/bam/logs/raw.erne-bs5.log 2>&1 || exit \$?

        samtools view -H unsorted.erne-bs5.bam > header.txt
        for ((i=1;i<=\$(grep -c SM:no_sample_specified <(samtools view -H unsorted.erne-bs5.bam));i++));
        do sed -i "0,/SM:no_sample_specified/{s/SM:no_sample_specified/SM:sample\${i}/}" header.txt; done || exit \$?

        samtools reheader header.txt unsorted.erne-bs5.bam |
        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o ${replicate}/bam/raw.erne-bs5.bam -
        """

}



// SENSITIVE READ ALIGNMENT USING "segemehl"
process "segemehl" {

    label 'finish'
    tag "$replicate"

    publishDir "${params.output}/bam", pattern: "$replicate/bam/raw.segemehl.bam", mode: 'copy', enabled: params.keepBams ? true : false
    publishDir "${params.output}/bam", pattern: "$replicate/bam/logs/raw.segemehl.log", mode: 'move'

    input:
    tuple val(replicate), val(readtype), path(reads)
    // eg. [replicate, ["input"], /path/to/inputs]
    path ctidx
    path gaidx
    path fasta
    path lamfa

    output:
    tuple val(replicate), path("$replicate/bam/raw.segemehl.bam")
    // eg. [replicate, /path/to/replicate]
    //path "$replicate/bam/raw.segemehl.bam"
    path "$replicate/bam/logs/raw.segemehl.log" 

    when:
    params.segemehl || ( params.merge && readtype == "merge" )

    script:
    if( params.SE )
        """
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs
        segemehl.x -i ${ctidx} -j ${gaidx} \\
        -d ${params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "${fasta}" : "${fasta} ${lamfa}"} \\
        -q ${reads} -o raw.segemehl.sam -I ${params.maxIns} -A ${params.minAccuracy} -s -t ${task.cpus} -F 1 -H 1 -D 1 \\
        > ${replicate}/bam/logs/raw.segemehl.log 2>&1 || exit \$?
        samtools view -Sb raw.segemehl.sam > ${replicate}/bam/raw.segemehl.bam
        """
    else
        """
        mkdir ${replicate} ${replicate}/bam ${replicate}/bam/logs
        segemehl.x -i ${ctidx} -j ${gaidx} \\
        -d ${params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? "${fasta}" : "${fasta} ${lamfa}"} \\
        -q ${reads[0]} -p ${reads[1]} -o raw.segemehl.sam -I ${params.maxIns} -A ${params.minAccuracy} -s -t ${task.cpus} -F 1 -H 1 -D 1 \\
        > ${replicate}/bam/logs/raw.segemehl.log 2>&1 || exit \$?
        samtools view -Sb raw.segemehl.sam > ${replicate}/bam/raw.segemehl.bam
        """

}



// POST-PROCESS BAMFILES READY FOR DOWNSTREAM ANALYSIS (erne-bs5)
process "erne_bs5_processing" {

    label 'finish'
    tag "$replicate"

    publishDir "${params.output}/bam", pattern: "$replicate/bam/proc.erne-bs5.bam", mode: 'copy', \
            enabled: params.keepBams || (!params.merge && params.noLambda && params.split == "${baseDir}/data/lambda.fa") ? true : false
    publishDir "${params.output}/bam", pattern: "${replicate}.bam", mode: 'copyNoFollow', \
            enabled: !params.merge && params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? true : false

    input:
    tuple val(replicate), path(erne)
    // eg. [replicate, /path/to/input]
    path fasta
    path lamfa
    
    output:
    tuple val(replicate), val("erne-bs5"), path("$replicate/bam/proc.erne-bs5.bam")
    // eg. [replicate, /path/to/replicate]
    path "${replicate}.bam"
    
    script:
    if (params.SE && params.noLambda && params.split == "${baseDir}/data/lambda.fa")
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s ${replicate}/bam/proc.erne-bs5.bam ${replicate}.bam

        samtools index ${erne}
        correct_sam_cigar.py ${erne} corrected.erne-bs5.bam || exit \$?
    
        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o sorted.erne-bs5.bam corrected.erne-bs5.bam
        samtools index sorted.erne-bs5.bam
        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . ${fasta} sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
    else if (!params.SE && params.noLambda && params.split == "${baseDir}/data/lambda.fa")
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s ${replicate}/bam/proc.erne-bs5.bam ${replicate}.bam

        samtools index ${erne}
        correct_sam_format.py -i ${params.maxIns} -t ${task.cpus} -T . ${erne} corrected.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -no unsorted.erne-bs5.bam corrected.erne-bs5.bam
        correct_sam_tlens -i unsorted.erne-bs5.bam > tlens.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o sorted.erne-bs5.bam tlens.erne-bs5.bam
        samtools index sorted.erne-bs5.bam
        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . ${fasta} sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
    else if (params.SE)
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s ${replicate}/bam/proc.erne-bs5.bam ${replicate}.bam

        if [ -z "\$(tail -c 1 ${fasta})" ]
        then
        cat ${fasta} ${lamfa} > fasta.tmp
        else
        cat ${fasta} > fasta.tmp; echo >> fasta.tmp; cat ${lamfa} >> fasta.tmp
        fi

        samtools index ${erne}
        correct_sam_cigar.py ${erne} corrected.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o sorted.erne-bs5.bam corrected.erne-bs5.bam
        samtools index sorted.erne-bs5.bam
        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . fasta.tmp sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
    else
        """
        mkdir ${replicate} ${replicate}/bam
        ln -s ${replicate}/bam/proc.erne-bs5.bam ${replicate}.bam

        if [ -z "\$(tail -c 1 ${fasta})" ]
        then
        cat ${fasta} ${lamfa} > fasta.tmp  
        else
        cat ${fasta} > fasta.tmp; echo >> fasta.tmp; cat ${lamfa} >> fasta.tmp      
        fi

        samtools index ${erne}
        correct_sam_format.py -i ${params.maxIns} -t ${task.cpus} -T . ${erne} corrected.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -no unsorted.erne-bs5.bam corrected.erne-bs5.bam
        correct_sam_tlens -i unsorted.erne-bs5.bam > tlens.erne-bs5.bam || exit \$?

        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o sorted.erne-bs5.bam tlens.erne-bs5.bam
        samtools index sorted.erne-bs5.bam

        filter_sam_erne.py -c ${params.XF} -t ${task.cpus} -T . fasta.tmp sorted.erne-bs5.bam ${replicate}/bam/proc.erne-bs5.bam
        """
}



// POST-PROCESS BAMFILES READY FOR DOWNSTREAM ANALYSIS (segemehl)
process "segemehl_processing" {

    label 'low'
    label 'finish'
    tag "$replicate"

    publishDir "${params.output}/bam", pattern: "$replicate/bam/proc.segemehl.bam", mode: 'copy', \
            enabled: params.keepBams || (!params.merge && params.noLambda && params.split == "${baseDir}/data/lambda.fa") ? true : false
    publishDir "${params.output}/bam", pattern: "${replicate}.bam", mode: 'copyNoFollow', \
            enabled: !params.merge && params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? true : false

    input:
    tuple val(replicate), path(sege)
    // eg. [replicate, /path/to/input]

    output:
    tuple val(replicate), val("segemehl"), path("$replicate/bam/proc.segemehl.bam")
    // eg. [replicate, segemehl, /path/to/proc.segemehl.bam]
    path "${replicate}.bam"

    script:
    if(params.SE)
        """
        mkdir ${replicate} ${replicate}/bam
        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o sorted.segemehl.bam ${sege}
        samtools index sorted.segemehl.bam

        filter_sam_xf_tag.py -c ${params.XF} sorted.segemehl.bam ${replicate}/bam/proc.segemehl.bam || exit \$?
        ln -s ${replicate}/bam/proc.segemehl.bam ${replicate}.bam
        """
    else
        """
        mkdir ${replicate} ${replicate}/bam
        correct_sam_tlens -i ${sege} > tlens.segemehl.bam || exit \$?
        samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
        -o sorted.segemehl.bam tlens.segemehl.bam
        samtools index sorted.segemehl.bam
        
        filter_sam_xf_tag.py -c ${params.XF} sorted.segemehl.bam ${replicate}/bam/proc.segemehl.bam || exit \$?
        ln -s ${replicate}/bam/proc.segemehl.bam ${replicate}.bam
        """
}



// MERGE BAMFILES READY FOR DOWNSTREAM ANALYSIS (OPTIONAL)
process "bam_merging" {

    label 'low'
    label 'finish'
    tag "$replicate"

    publishDir "${params.output}/bam", pattern: "$replicate/bam/merged.bam", mode: 'copy', \
            enabled: params.keepBams || (params.noLambda && params.split == "${baseDir}/data/lambda.fa") ? true : false
    publishDir "${params.output}/bam", pattern: "${replicate}.bam", mode: 'copyNoFollow', \
            enabled: params.noLambda && params.split == "${baseDir}/data/lambda.fa" ? true : false

    input:
    tuple val(replicate), val(erne), path(erne_bs5), val(sege), path(segemehl)
    // eg. [replicate, proc.erne-bs5.bam, proc.segemehl.bam]

    output:
    tuple val(replicate), val("merged"), path("$replicate/bam/merged.bam")
    // eg. [replicate, merged, /path/to/replicate/replicate.bam]
    path "${replicate}.bam"
   
    when:
    params.merge
    
    script:
    """
    mkdir ${replicate} ${replicate}/bam
    ln -s ${replicate}/bam/merged.bam ${replicate}.bam
    samtools merge ${replicate}/bam/merged.bam ${erne_bs5} ${segemehl}
    """
}



// SPLIT BAMFILES FOR LAMBDA AND EVERTHING ELSE
process "bam_subsetting" {

    label 'low'
    label 'finish'
    tag "$replicate"

    publishDir "${params.output}/bam", pattern: "${replicate}/bam/${chrom}.bam", mode: 'copy', enabled: params.keepBams ? true : false
    publishDir "${params.output}/bam", pattern: "${replicate}/bam/subset.bam", mode: 'copy', enabled: true
    publishDir "${params.output}/bam", pattern: "${replicate}.bam", mode: 'copyNoFollow', enabled: true

    input:
    tuple val(replicate), val(bamtype), path(bamfile)
    // eg. [replicate, merged, bamfile.bam]
    path fai
    path lai
    val chrom

    output:
    tuple val(replicate), val("lambda"), path("${replicate}/bam/${chrom}.bam")
    tuple val(replicate), val("subset"), path("${replicate}/bam/subset.bam")
    // eg. [replicate, subset, /path/to/replicate/bam/subset.bam]
    path "${replicate}.bam"
    
    when:
    !params.noLambda || params.split != "${baseDir}/data/lambda.fa"
   
    script:
    """
    mkdir ${replicate} ${replicate}/bam
    ln -s ${replicate}/bam/subset.bam ${replicate}.bam

    # reset alignments whose mates have erroneously aligned to lambda
    samtools view -h ${bamfile} |
    awk 'BEGIN {OFS="\\t"} {if((\$3=="${chrom}" && \$7=="=") || \$7=="${chrom}") \
    {if(\$2==97){\$2=73} else { if(\$2==81){\$2=89} \
    else { if(\$2==161){\$2=137} else{ if(\$2==145){\$2=153} }}}; \$7="*"; \$8=0}; {print \$0}}' |
    
    samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
    -o sort.bam - || exit \$?
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

    publishDir "${params.output}/bam", pattern: "${replicate}/${replicate}.bam.stats", mode: 'copy'
    publishDir "${params.output}/bam", pattern: "${replicate}/stats/*.png", mode: 'move'

    input:
    tuple val(replicate), val(bamtype), path(bamfile)
    // eg. [replicate, bamtype, bamfile.bam]

    output:
    path "${replicate}/${replicate}.bam.stats"
    path "${replicate}/stats/*.png"    

    script:
    """
    mkdir ${replicate} ${replicate}/stats
    samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
    -o sorted.bam ${bamfile} || exit \$?
    samtools stats sorted.bam > ${replicate}/${replicate}.bam.stats || exit \$?
    plot-bamstats -p ${replicate}/stats/ ${replicate}/${replicate}.bam.stats
    """    
}



// PRE-PROCESS BAMFILES READY FOR DOWNSTREAM METHYLATION EXTRACTION
process "bam_filtering" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"

    publishDir "${params.output}/bam", pattern: "$replicate/bam/*.bam", mode: 'copy', enabled: params.keepBams && "$bamtype" != "lambda" ? true : false
    publishDir "${params.output}/bam", pattern: "${replicate}.bam", mode: 'copyNoFollow', enabled: true

    input:
    tuple val(replicate), val(bamtype), path(bamfile)
    // eg. [replicate, lambda, /path/to/bamfile.bam]
    // eg. [replicate, subset, /path/to/bamfile.bam]

    output:
    tuple val(replicate), val(bamtype), path("$replicate/bam/unique.bam")
    tuple val(replicate), val(bamtype), path("$replicate/bam/*.bam")
    tuple val(replicate), val(bamtype), path("${replicate}.bam")
    // eg. [replicate, lambda, /path/to/replicate/*.bam]
    // eg. [replicate, subset, /path/to/replicate/*.bam]
    
    when:
    params.unique

    script:
    """
    mkdir ${replicate} ${replicate}/bam
    ln -s ${replicate}/bam/unique.bam ${replicate}.bam
    filter_sam_uniqs.py ${bamfile} ${replicate}/bam/unique.bam ${replicate}/bam/multimapped.bam
    """
}
