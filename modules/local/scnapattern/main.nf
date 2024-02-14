

process CALCULATE_SCNAPATTERN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "dincalcilab/pandas-pyranges:0.0.111-f05923d"

    input:
    tuple val(meta), path(segmentfile)

    output:
    tuple val(meta), path("*_classification.txt"), emit: classification_table
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = "${meta.format}"
    def ploidy = "${meta.ploidy}"

    """

    calculate_scnapattern.py \\
        --ploidy $ploidy \\
        --file-format $format \\
        --sample-name $prefix \\
        $args
        $segmentfile
        ${prefix}_classification.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scnapattern: \$(calculate_scnapattern.py --version |& sed '1!d ; s/calculate_scnapattern.py //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = "${meta.format}"
    def ploidy = "${meta.ploidy}"

    """

    touch ${prefix}_classification.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scnapattern: \$(calculate_scnapattern.py --version |& sed '1!d ; s/calculate_scnapattern.py //')
    END_VERSIONS
    """
}
