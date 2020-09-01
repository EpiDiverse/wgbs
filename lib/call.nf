// SPLIT BAMS ACCORDING TO READ GROUPS
process "bam_grouping" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"

    //publishDir "${params.output}/bam"

    input:
    tuple val(replicate), val(bamtype), path(bam)
    // eg. [replicate, lambda, /path/to/unsorted.bam]
    // eg. [replicate, subset, /path/to/unsorted.bam]

    output:
    tuple val(replicate), path("*.txt")
    // eg. [replicate, no_sample_specified.txt]
    // eg. [replicate, [sample1.txt, sample2.txt, ... sampleN.txt]]
    
    when:
    !params.ignoreRG && bamytype != "lambda"

    script:
    """
    samtools view -H ${bam} | grep "^@RG" |
    awk '{for(i=1;i<=NF;i++){if(\$i~"^ID"){split(\$i,ID,":")}else{if(\$i~"^SM"){split(\$i,SM,":")}}};
    print ID[2] >> SM[2]".txt"}'
    """ 
}

// [replicate, [sample1.txt, sample2.txt, ... sampleN.txt]]
// txts.transpose().map{ tuple(it[0], it[1].tokenize(".").init().join(""), it[1]) }
// [replicate, sample1, sample1.txt], [replicate, sample2, sample2.txt], ...

// [replicate, subset, /path/to/unsorted.bam]
// bams.filter( it[1] == "subset" ).combine(txts, by:0)
// [replicate, subset, /path/to/unsorted.bam, sample1, sample1.txt]

// PRE-PROCESS BAMFILES READY FOR DOWNSTREAM METHYLATION EXTRACTION
process "bam_sampling" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"

    //publishDir "${params.output}/bam"

    input:
    tuple val(replicate), val(bamtype), path("input.bam"), val(filename), path("input.txt")
    // eg. [replicate, subset, /path/to/unsorted.bam, sample1, sample1.txt]

    output:
    tuple val(replicate), val(bamtype), path("sample.bam"), val(filename)
    // eg. [replicate, subset, /path/to/sample.bam, sample1]
    
    when:
    !params.ignoreRG && bamtype != "lambda"

    script:
    """
    samtools view -bhR input.txt input.bam > sample.bam
    """ 
}


// bams.filter{ it[1] == "lambda" }.map{ tuple(it[0], it[1], it[0], it[2]) }

// PRE-PROCESS BAMFILES READY FOR DOWNSTREAM METHYLATION EXTRACTION
process "bam_processing" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"

    //publishDir "${params.output}/bam"

    input:
    tuple val(replicate), val(bamtype), path("unsorted.bam"), val(filename)
    // eg. [replicate, lambda, /path/to/unsorted.bam, replicate]
    // eg. [replicate, subset, /path/to/unsorted.bam, sample1]

    output:
    tuple val(replicate), val(bamtype), val(filename), path("unique.bam")
    // eg. [replicate, lambda, replicate, /path/to/unique.bam]
    // eg. [replicate, subset, sample1, /path/to/unique.bam]
    
    script:
    """
    samtools sort -T deleteme -o sorted.bam unsorted.bam
    change_sam_qname -i sorted.bam -o unique.bam --tags ${params.SE ? "HI" : "HI XB"} --read_name_tag XN
    """ 
}


// MARK DUPLICATES (OPTIONAL)
process "Picard_MarkDuplicates" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"

    publishDir "${params.output}/bam", pattern: "$replicate/bam/*.bam", mode: 'copy', enabled: params.keepBams ? true : false
    publishDir "${params.output}/bam", pattern: "$replicate/duplicates*.txt", mode: 'move'
    publishDir "${params.output}/bam", pattern: "$replicate/*/logs/*.log", mode: 'move'

    input:
    tuple val(replicate), val(bamtype), val(filename), path(bam)
    // eg. [replicate, lambda, replicate, /path/to/*.bam] or [replicate, subset, sample1, /path/to/*.bam]

    output:
    tuple val(replicate), val(bamtype), val(filename), path("$replicate/*/*.bam")
    tuple val(replicate), val(bamtype), val(filename), path("$replicate/*.txt")
    path "$replicate/*/logs/*.log"

    when:
    !params.noDedup

    script:
    duplicates = "${bamtype == "lambda" ? "lambda" : "duplicates"}"
    """
    mkdir tmp ${replicate} ${replicate}/lambda ${replicate}/bam ${replicate}/bam/logs
    
    picard -Xmx${task.memory.getBytes() - 2147483648} MarkDuplicates TMP_DIR=tmp \\
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=\$(ulimit -n) \\
    VALIDATION_STRINGENCY=LENIENT \\
    I=${bam} O=${replicate}/${bamtype == "lambda" ? "lambda" : "bam"}/${replicate == filename ? "markDups" : "markDups.${filename}"}.bam \\
    M=${replicate}/${replicate == filename ? "${duplicates}" : "${duplicates}.${filename}"}.txt \\
    > ${replicate}/bam/logs/markDups.${replicate == filename || bamtype == "lambda" ? "${bamtype}" : "${filename}"}.log 2>&1
    """
}



// METHYLATION CALLING USING "MethylDackel"
process "MethylDackel" {

    label 'low'
    label 'ignore'
    tag "$replicate - $bamtype"

    publishDir "${params.output}/bedGraph", pattern: "CpG/*.bedGraph", mode: 'copy'
    publishDir "${params.output}/bedGraph", pattern: "CHG/*.bedGraph", mode: 'copy'
    publishDir "${params.output}/bedGraph", pattern: "CHH/*.bedGraph", mode: 'copy'
    publishDir "${params.output}/bam", pattern: "$replicate/*.svg", mode: 'move'
    publishDir "${params.output}/bedGraph", pattern: "logs/*.err", mode: 'move'

    input:
    tuple val(replicate), val(bamtype), val(filename), path(bam)
    // eg. [replicate, lambda, replicate, markDups.bam]
    // eg. [replicate, subset, sample1, markDups.bam]
    path fasta
    path lamfa
    val context
    
    output:
    tuple val(replicate), val(bamtype), val(filename), path("*/*.bedGraph")
    tuple val(replicate), val(bamtype), val(filename), path("*/*.svg")
    path "logs/*.err"

    script:
    """
    mkdir logs bedGraph ${replicate} ${bamtype}
    samtools index ${bam}

    STR=\$(echo \$(MethylDackel mbias ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} ${bam} \\
    ${bamtype == "lambda" ? "lambda" : "$replicate"}/${filename} ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"} 2>&1 | cut -d ":" -f2))
    MethylDackel extract ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \\
    ${bam} ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"}-o ${bamtype == "lambda" ? "lambda" : "$replicate"}/${filename} \$STR \\
    > logs/${bamtype}.${replicate == filename ? "${replicate}" : "${replicate}.${filename}"}.err 2>&1

    find ${replicate} -name "*.bedGraph" -type f |
    while read file; do id=\$(basename \$file .bedGraph);
    mkdir \${id##*_}; mv \$file \${id##*_}/${replicate == filename ? "${replicate}" : "${replicate}.${filename}"}.bedGraph;
    done
    """
}


// LINEAR REGRESSION OF DUPLICATES
process "linear_regression" {

    label 'low'
    label 'ignore'

    publishDir "${params.output}", pattern: "Duplicates.*", mode: 'move'

    input:
    path "duplicates"
    // eg. [/path/to/duplicates.txt, ...]

    output:
    path "Duplicates.*"

    when:
    !params.noDedup

    script:
    """
    ls duplicates* | while read file; do grep -A2 "^## METRICS CLASS" \$file | 
    tail -1 | cut -f ${params.SE ? "2,6" : "3,7"} >> Duplicates.tsv; done

    Rscript ${baseDir}/bin/scatterPlot.R Duplicates.tsv
    """
}


// ESTIMATION OF CONVERSION RATE FROM LAMBDA
process "conversion_rate_estimation" {

    label 'low'
    label 'ignore'
    tag "$replicate - ${ bamtype == "lambda" ? "${chrom}" : "${params.chrom}" }"
    
    input:
    tuple val(replicate), val(bamtype), val(filename), path("bedGraph")
    // eg. [replicate, lambda, replicate, [CpG.bedGraph, CHG.bedGraph, CHH.bedGraph]]
    // eg. [replicate, subset, sample1, [CpG.bedGraph, CHG.bedGraph, CHH.bedGraph]]
    // lambda bamtype will always have each context in bedGraph files
    val chrom

    output:
    tuple val(replicate), path("*.txt")
       
    when:
    ( bamtype == "lambda" || params.chrom )

    script:
    """
    echo -e "${replicate}\\t${bamtype == "lambda" ? "${chrom}" : "${params.chrom}"}\\tNon-conversion Rate (%): \\
    \$(tail -q -n+2 bedGraph* | awk '\$1~/${bamtype == "lambda" ? "${chrom}" : "${params.chrom}"}/{{m += \$5}; {u += \$6}} END {t = (m+u); print (m/t)*100}')" \\
    > ${replicate}.txt
    """ 
}










/*
// PRE-PROCESS BAMFILES READY FOR DOWNSTREAM METHYLATION EXTRACTION
process "bam_processing" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"
    
    input:
    tuple replicate, bamtype, path("unsorted.bam")
    // eg. [replicate, lambda, /path/to/unsorted.bam]
    // eg. [replicate, subset, /path/to/unsorted.bam]

    output:
    tuple replicate, bamtype, path("unique.bam")
    // eg. [replicate, lambda, /path/to/unique.bam]
    // eg. [replicate, subset, /path/to/unique.bam]
    
    script:
    """
    samtools sort -T deleteme -o sorted.bam unsorted.bam
    change_sam_qname -i sorted.bam -o unique.bam --tags ${params.SE ? "HI" : "HI XB"} \\
    --read_name_tag XN
    """ 
}


// MARK DUPLICATES (OPTIONAL)
process "Picard_MarkDuplicates" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"

    input:
    tuple replicate, bamtype, path(bam)
    // eg. [replicate, lambda, /path/to/*.bam] or [replicate, subset, /path/to/*.bam]

    output:
    tuple replicate, bamtype, path("$replicate/bam/markDups.bam")
    tuple replicate, bamtype, path("$replicate/duplicates.txt")
    path "$replicate/bam/logs/*.log"

    when:
    !params.noDedup

    script:
    """
    mkdir tmp ${replicate} ${replicate}/bam ${replicate}/bam/logs
    
    picard -Xmx${task.memory.getBytes() - 2147483648} MarkDuplicates TMP_DIR=tmp \\
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=\$(ulimit -n) \\
    VALIDATION_STRINGENCY=LENIENT \\
    I=${bam} O=${replicate}/bam/markDups.bam M=${replicate}/duplicates.txt \\
    > ${replicate}/bam/logs/markDups.${bamtype}.log 2>&1
    """
}



// METHYLATION CALLING USING "MethylDackel"
process "MethylDackel" {

    label 'low'
    label 'ignore'
    tag "$replicate - $bamtype"

    input:
    tuple replicate, bamtype, path(bam)
    // eg. [replicate, lambda, markDups.bam]
    path fasta
    path lamfa
    val context
    
    output:
    tuple replicate, bamtype, path("* /*.bedGraph")
    tuple replicate, bamtype, path("$replicate/*.svg")
    path "logs/*.err"

    script:
    """
    mkdir logs bedGraph ${replicate}
    BAM=\$(ls *.bam)
    samtools index \$BAM

    STR=\$(echo \$(MethylDackel mbias ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \$BAM ${replicate}/Mbias ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"} 2>&1 | cut -d ":" -f2))
    MethylDackel extract ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \\
    \$BAM ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"}-o ${replicate}/${replicate} \$STR \\
    > logs/${bamtype}.${replicate}.err 2>&1

    find ${replicate} -name "*.bedGraph" -type f |
    while read file; do id=\$(basename \$file .bedGraph);
    mkdir \${id##*_}; mv \$file \${id##*_}/${replicate}.bedGraph;
    done
    """
}


// LINEAR REGRESSION OF DUPLICATES
process "linear_regression" {

    label 'low'
    label 'ignore'

    input:
    path "duplicates"
    // eg. [/path/to/duplicates.txt, ...]

    output:
    path "Duplicates.*"

    when:
    !params.noDedup

    script:
    """
    ls duplicates* | while read file; do grep -A2 "^## METRICS CLASS" \$file | 
    tail -1 | cut -f ${params.SE ? "2,6" : "3,7"} >> Duplicates.tsv; done

    Rscript ${baseDir}/bin/scatterPlot.R Duplicates.tsv
    """
}


// ESTIMATION OF CONVERSION RATE FROM LAMBDA
process "conversion_rate_estimation" {

    label 'low'
    label 'ignore'
    tag "$replicate - ${ bamtype == "lambda" ? "${chrom}" : "${params.chrom}" }"
    
    input:
    tuple replicate, bamtype, path("bedGraph")
    // eg. [replicate, lambda, [CpG.bedGraph, CHG.bedGraph, CHH.bedGraph]]
    // eg. [replicate, subset, [CpG.bedGraph, CHG.bedGraph, CHH.bedGraph]]
    // lambda bamtype will always have each context in bedGraph files
    val chrom

    output:
    tuple replicate, path("*.txt")
       
    when:
    ( bamtype == "lambda" || params.chrom )

    script:
    """
    echo -e "${replicate}\\t${bamtype == "lambda" ? "${chrom}" : "${params.chrom}"}\\tNon-conversion Rate (%): \\
    \$(tail -q -n+2 bedGraph* | awk '\$1~/${bamtype == "lambda" ? "${chrom}" : "${params.chrom}"}/{{m += \$5}; {u += \$6}} END {t = (m+u); print (m/t)*100}')" \\
    > ${replicate}.txt
    """ 
}
*/