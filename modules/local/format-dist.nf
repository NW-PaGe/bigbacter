process FORMAT_DIST {
    tag "${taxa}-${cluster}"
    label 'process_low'

    input:
    tuple val(taxa), val(cluster), path(dist), path(tree)

    output:
    tuple val(taxa), val(cluster ), path("dist.formatted.txt"), path(tree), emit: dist

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    """    
    # filter dist file to contain only samples in the core SNP tree.
    # doing this with shell is much faster than R or Python
    format-dist.sh "${dist}" "${tree}" ${dist.name.endsWith('.csv.gz') ? 'bb' : 'pp'}
    """
}