process BB_CLUSTER {
    tag "${taxa}"
    label 'process_high'

    input:
    tuple val(taxa), path(assembly), path(db)
    val timestamp

    output:
    tuple val(taxa), path('clusters.csv'),    emit: clusters    
    tuple val(taxa), path('*.tar.gz'),        emit: db
    tuple val(taxa), path("dist.all.csv.gz"), emit: dist, optional: true
    tuple val(taxa), path("mds.dot"),         emit: dot,  optional: true
    tuple val(taxa), path("tree.nwk"),        emit: tree, optional: true
    path 'versions.yml',                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    prefix = "${timestamp}-${taxa}"
    """
    # Prepare current database files
    if [ -e "${db}" ]
    then
        ## Extract current db contents
        tar -xzvf ${db} -C ./
        ## Set paths for signature files and cluster info
        DB_FILES="--sigs */sig/*.sig --clusters */clusters.csv"
    else
        DB_FILES=''
    fi

    # Assign clusters
    bb-cluster.py \\
        ${args} \\
        --dist ${params.bb_dist} \\
        --outdir new_db \\
        --query ${assembly} \\
        \${DB_FILES}

    # Format outputs
    ## Prepare new database
    mkdir ${timestamp}
    mv new_db/sig/ ${timestamp}/
    cp new_db/clusters.csv ${timestamp}/
    tar -czvf ${timestamp}.tar.gz ${timestamp}

    ## Prepare files for publishing
    mv new_db/* ./ || true
    gzip dist.all.csv || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bb-cluster.py: \$(bb-cluster --version )
    END_VERSIONS
    """
}
