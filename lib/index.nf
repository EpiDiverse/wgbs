#!/usr/bin/env nextflow


// FAST READ ALIGNMENT USING "erne-bs5" 
process "erne_bs5_indexing" {

    label "low"
    label "finish"

    publishDir "${params.output}/index", mode: 'copy', enabled: params.INDEX ? true : false

    input:
    path fasta
    path fai
    path lamfa
    path lai

    output:
    path "*.ebm"

    when:
    !params.segemehl || params.merge

    script:
    if( params.noLambda && params.split == "${baseDir}/data/lambda.fa" )
        """
        erne-create --methyl-hash --output-prefix ${fasta.baseName} --fasta ${fasta}
        """
    else
        """
        erne-create --methyl-hash --output-prefix ${fasta.baseName} --fasta ${fasta} --fasta ${lamfa}
        """
}


// INDEXING FOR "erne-bs5" 
process "segemehl_indexing" {

    label "low"
    label "finish"

    publishDir "${params.output}/index", mode: 'copy', enabled: params.INDEX ? true : false

    input:
    path fasta
    path fai
    path lamfa
    path lai

    output:
    path "*.ctidx"
    path "*.gaidx"

    when:
    params.segemehl || params.merge

    script:
    if( params.noLambda && params.split == "${baseDir}/data/lambda.fa" )
        """
        segemehl.x -x ${fasta.baseName}.ctidx -y ${fasta.baseName}.gaidx -d ${fasta} -F 1
        """
    else
        """
        segemehl.x -x ${fasta.baseName}.ctidx -y ${fasta.baseName}.gaidx -d ${fasta} ${lamfa} -F 1
        """
}