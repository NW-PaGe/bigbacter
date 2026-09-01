process DS_RATE {

    input:
    tuple val(meta), path(jsons), path(ref)

    output:
    tuple val(meta), path("ds_rate.csv"), emit: result
    path 'versions.yml',                  emit: versions

    script:
    tool = 'bigbacter_ds_rate.py'
    """
    ${tool} \\
        --ref ${ref} \\
        --max-depth ${params.max_depth} \\
        ${jsons}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}