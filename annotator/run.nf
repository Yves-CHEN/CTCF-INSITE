#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.chipF=""
params.genomeRef=""


process runRscript {
    container = 'docker://yveschen/ctcfinsite-annotator:v1.0'
    cpus 1
    publishDir 'results', mode: 'copy', pattern: 'results/*.fullAnnot.txt'

    input:
        tuple  path(chipFile), path(functionFile),path(rscript), path(resource), path(ref)
     output:
        path "results/*.fullAnnot.txt"
    script:
    """
    Rscript $rscript $chipFile  ${task.cpus} ${ref}
    """


}

workflow {
    paras = tuple(file(params.chipF), file("$baseDir/common.R"), file("$baseDir/annotate.R"), file ("$baseDir/resource"), file(params.genomeRef))
    runRscript(paras)
}

