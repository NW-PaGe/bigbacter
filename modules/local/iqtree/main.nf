process IQTREE {
    label 'process_high'
    
    input:
    tuple val(meta), path(aln), path(const_sites)

    output:
    tuple val(meta), path(treefile), emit: result
    path 'versions.yml',             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args     = task.ext.args ?: ''
    prefix   = "${meta.timestamp}-${meta.taxa}-${meta.cluster}"
    bs       = "${meta.n > 4 ? '-B 1000' : ''}"
    treefile = "${meta.n > 4 ? '*.contree' : '*.treefile'}"
    """
    # run IQTREE2
    iqtree3 \\
        -s ${aln} \\
        -fconst \$(cat ${const_sites}) \\
        -T ${task.cpus} \\
        ${args} \\
        ${bs}
    
    #### VERSION INFO ####
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree3 --version 2>&1 | head -n 1 | sed 's/^IQ-TREE multicore version //; s/ .*//')
    END_VERSIONS
    """
}