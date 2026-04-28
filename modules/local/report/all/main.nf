process ALL {
    label 'process_low'
        
    input:
    tuple val(meta), path(summaries, stageAs: "input/*")

    output:
    tuple val(meta), path("*.csv"), emit: result
    path "versions.yml", emit: versions

    script:
    tool = 'bigbacter_report_all.py'
    """
    ${tool} \\
        --timestamp "${meta.timestamp}" \\
        input/


    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
