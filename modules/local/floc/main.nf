process FLOC {
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(meta), path(assemblies, stageAs: 'input/*'), path(db, stageAs: 'input/')

    output:
    tuple val(meta), path("clusters.csv"), emit: clusters
    tuple val(meta), path("sigs/*"),       emit: sigs
    tuple val(meta), path("dist/*"),       emit: dist
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' 
    tool = 'floc'
    """
    ${tool} \\
        --pairwise \\
        -d ${params.clust_dist} \\
        --abund-z ${params.clust_z} \\
        --min-global-abund ${params.clust_min_global} \\
        --max-global-dist ${params.clust_max_global} \\
        ${params.clust_ignore_qc ? '--ignore-qc' : ''} \\
        ${args} \\
        input/

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}