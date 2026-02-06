process CLUSTER {
    label 'process_low'
        
    input:
    tuple val(meta), path(input, stageAs: "input/*")
    path(samplesheet, stageAs: 'input/samplesheet.csv')
    path(microreact_template, stageAs: 'input/template.microreact')

    output:
    tuple val(meta), path("*.microreact"), emit: report
    tuple val(meta), path("*_summary.*"),  emit: summary
    tuple val(meta), path("*-dist.*"),     emit: dist
    tuple val(meta), path("*.nwk"),        emit: tree, optional: true
    path "versions.yml", emit: versions

    script:
    tool = 'bigbacter_report.py'
    """
    ${tool} \\
        --timestamp "${meta.timestamp}" \\
        --taxa "${meta.taxa}" \\
        --cluster "${meta.cluster}" \\
        --strong-link ${params.strong_link_threshold} \\
        --inter-link ${params.inter_link_threshold} \\
        --partition-distance ${params.partition_distance} \\
        ${ meta.recomb_masked ? "--recomb-masked" : '' } \\
        input/


    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
