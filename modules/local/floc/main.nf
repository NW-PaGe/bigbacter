process FLOC {
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(meta), path(assemblies, stageAs: 'input/asm/*'), path(db, stageAs: 'input/')

    output:
    tuple val(meta), path("clusters.csv"),           emit: clusters
    tuple val(meta), path("sigs/*"),                 emit: sigs
    tuple val(meta), path("dist/*"),                 emit: dist
    tuple val(meta), path("global_containment.csv"), emit: qc
    tuple val(meta), path("nj_tree.html"),           emit: plot, optional: true
    tuple val(meta), path("nj_tree.nwk"),            emit: tree, optional: true
    tuple val(meta), path("config.json"),            emit: config, optional: true


    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' 
    tool = 'floc'
    // config.json is only written when the taxa database is created (no existing sig db)
    def new_db = !db
    """
    ${tool} \\
        --pairwise \\
        -d ${params.clust_dist} \\
        --min-hash-frac ${params.clust_min_hash_frac} \\
        --min-hash-freq ${params.clust_min_hash_freq} \\
        --assign-hash-freq ${params.clust_assign_hashfreq} \\
        ${params.clust_ignore_qc ? '--ignore-qc' : ''} \\
        ${params.clust_overwrite ? '--overwrite' : ''} \\
        ${params.clust_plot ? '--plot' : ''} \\
        ${args} \\
        input/*/

    # record the settings a new taxa database is created with (pushed to the database with --push).
    # existing databases already have a config.json, so none is written.
    if ${new_db}; then
        cat <<-END_CONFIG > config.json
    {
      "clust_dist": ${params.clust_dist},
      "clust_assign_hashfreq": ${params.clust_assign_hashfreq},
      "floc_version": "\$(${tool} --version 2>&1 | tr -d '\\r')",
      "bigbacter_version": "${workflow.manifest.version}"
    }
    END_CONFIG
    fi

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}