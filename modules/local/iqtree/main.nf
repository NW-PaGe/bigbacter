process IQTREE {
    label 'process_high'

    input:
    tuple val(meta), path(aln), path(const_sites)

    output:
    tuple val(meta), path("*.nwk")    , emit: result
    tuple val(meta), path("*.treefile"), emit: mltree
    tuple val(meta), path("*.contree"), optional: true, emit: contree
    tuple val(meta), path("*.iqtree") , emit: report
    tuple val(meta), path("*.log")    , emit: log
    path 'versions.yml'               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.timestamp}-${meta.taxa}-${meta.cluster}"
    def n      = (meta.n ?: 0) as Integer
    def bs     = n >= 4 ? '-B 1000' : ''
    """
    set -euo pipefail

    # -fconst expects a comma-separated list on the command line, not a file.
    # Strip all whitespace/newlines so the shell passes a single token.
    FCONST=\$(tr -d '[:space:]' < ${const_sites})
    if [ -z "\${FCONST}" ]; then
        echo "ERROR: ${const_sites} is empty" >&2
        exit 1
    fi

    iqtree3 \\
        -s ${aln} \\
        -fconst "\${FCONST}" \\
        -T ${task.cpus} \\
        --prefix ${prefix} \\
        --redo \\
        ${args} \\
        ${bs}

    # With -B, IQ-TREE writes both .treefile (ML) and .contree (UFBoot consensus).
    # Emit one deterministic file so downstream processes don't branch on n.
    if [ -s "${prefix}.contree" ]; then
        cp "${prefix}.contree" "${prefix}.nwk"
    else
        cp "${prefix}.treefile" "${prefix}.nwk"
    fi

    #### VERSION INFO ####
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree3 --version 2>&1 | head -n 1 | sed 's/^IQ-TREE multicore version //; s/ .*//')
    END_VERSIONS
    """
}