process SELECT_REF {

    input:
    tuple val(meta), path(assemblies, stageAs: 'input/*')

    output:
    tuple val(meta), path("ref.fa.gz"), path("ref.json"), emit: result
    path 'versions.yml',                                  emit: versions

    script:
    tool = 'bigbacter_select_ref.py'
    """
    ${tool} \\
        ${assemblies} \\
        --min-contig-len ${params.ref_min_contig_len} \\
        --contig-penalty ${params.ref_contig_penalty} \\
        --ksize ${params.ref_ksize} \\
        --scaled ${params.ref_scaled}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}