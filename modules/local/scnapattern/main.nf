

process SCNAPATTERN_MAIN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scnapattern: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scnapattern: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
