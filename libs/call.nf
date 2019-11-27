// PRE-PROCESS BAMFILES READY FOR DOWNSTREAM METHYLATION EXTRACTION
process "bam_sorting" {

    label 'low'
    label 'finish'
    tag "$replicate - $bamtype"
    
    input:
    tuple replicate, bamtype, path(bamfile)
    // eg. [replicate, lambda, /path/to/bamfile.bam]
    // eg. [replicate, subset, /path/to/bamfile.bam]

    output:
    tuple replicate, bamtype, path("$replicate/bam/*.bam")
    // eg. [replicate, lambda, /path/to/replicate/*.bam]
    // eg. [replicate, subset, /path/to/replicate/*.bam]
    
    script:
    """
    mkdir ${replicate} ${replicate}/bam
    samtools sort -T deleteme -o ${replicate}/bam/sorted.bam ${bamfile}
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
    tuple replicate, bamtype, path("$replicate/stats/duplicates.txt")
    path "$replicate/bam/logs/*.log"

    when:
    !params.noDedup

    script:
    """
    mkdir tmp ${replicate} ${replicate}/stats ${replicate}/bam ${replicate}/bam/logs
    picard -Xmx${task.memory.getBytes() - 2147483648} MarkDuplicates TMP_DIR=tmp \\
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=\$(ulimit -n) \\
    VALIDATION_STRINGENCY=LENIENT \\
    I=${bam} O=marked.bam M=${replicate}/stats/duplicates.txt \\
    > ${replicate}/bam/logs/markDups.${bamtype}.log 2>&1 || exit \$?
    change_sam_qname -i marked.bam -o ${replicate}/bam/markDups.bam --tags HI XB --read_name_tag XN --restore
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
    tail -1 | cut -f ${params.SE ? "2,6" : "3,7"} >> input.tsv; done

    Rscript ${baseDir}/bin/linearModel.R input.tsv
    """
}


// METHYLATION CALLING USING "MethylDackel"
process "MethylDackel" {

    label 'low'
    label 'ignore'
    tag "$replicate - $bamtype"

    input:
    tuple replicate, bamtype, path(bamfile)
    // eg. [replicate, lambda, markDups.bam]
    path fasta
    path lamfa
    val context
    
    output:
    tuple replicate, bamtype, path("$replicate/bedGraph/*.bedGraph")
    tuple replicate, bamtype, path("$replicate/stats/*.svg")
    path "$replicate/bedGraph/logs/*.err"

    script:
    if( !params.unique && !params.noDedup && ( params.segemehl || params.merge ))
        """
        mkdir ${replicate} ${replicate}/stats ${replicate}/bedGraph ${replicate}/bedGraph/logs
        BAM=\$(ls *.bam)
        change_sam_qname -i \$BAM -o restored.bam --tags HI XB --read_name_tag XN || exit \$?
        samtools index restored.bam

        STR=\$(echo \$(MethylDackel mbias ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} restored.bam ${replicate}/stats/Mbias ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"} 2>&1 | cut -d ":" -f2))
        MethylDackel extract ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \\
        restored.bam ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"}-o ${replicate}/bedGraph/${replicate} \$STR \\
        > ${replicate}/bedGraph/logs/${bamtype}.${replicate}.err 2>&1
        """
    else
        """
        mkdir ${replicate} ${replicate}/stats ${replicate}/bedGraph ${replicate}/bedGraph/logs
        BAM=\$(ls *.bam)
        samtools index \$BAM

        STR=\$(echo \$(MethylDackel mbias ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \$BAM ${replicate}/stats/Mbias ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"} 2>&1 | cut -d ":" -f2))
        MethylDackel extract ${bamtype == "lambda" ? "${lamfa}" : "${fasta}"} \\
        \$BAM ${bamtype == "lambda" ? "--CHH --CHG " : "${context}"}-o ${replicate}/bedGraph/${replicate} \$STR \\
        > ${replicate}/bedGraph/logs/${bamtype}.${replicate}.err 2>&1
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