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
    change_sam_qname -i sorted.bam -o unique.bam --tags HI ${params.SE ? "" : "XB"}\\
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
    tuple replicate, bamtype, path("bedGraph/*/*.bedGraph")
    tuple replicate, bamtype, path("$replicate/*.svg")
    path "bedGraph/logs/*.err"

    script:
    """
    mkdir logs ${replicate}
    BAM=\$(ls *.bam)
    samtools index \$BAM

    STR=\$(echo \$(MethylDackel mbias ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \$BAM ${replicate}/Mbias ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"} 2>&1 | cut -d ":" -f2))
    MethylDackel extract ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \\
    \$BAM ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"}-o ${replicate}/${replicate} \$STR \\
    > bedGraph/logs/${bamtype}.${replicate}.err 2>&1

    find ${replicate} -name "*.bedGraph" -type f |
    while read file; do id=\$(basename \$file .bedGraph);
    mkdir bedGraph/\${id##*_}; mv \$file bedGraph/\${id##*_}/${replicate}.bedGraph;
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
    tuple replicate, bamtype, path(bedGraph)
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
    \$(tail -q -n+2 ${bedGraph} | awk '\$1~/${bamtype == "lambda" ? "${chrom}" : "${params.chrom}"}/{{m += \$5}; {u += \$6}} END {t = (m+u); print (m/t)*100}')" \\
    > ${replicate}.txt
    """ 
}