process PUSH_CLUSTER_FILES {
    tag "${taxa}_${cluster}"
    label 'process_low'

    input:
    tuple val(taxa), val(cluster), path(ref, stageAs: "ref.fa.gz"), path(new_snippy), path(assembly)

    output:
    tuple path("*.tar.gz", includeInputs: true), path("*_T*.fa.gz", includeInputs: true), path('ref.fa.gz', includeInputs: true), emit: cluster_files
    
    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
    # apply compression if needed
    gzip * || true
    '''
}

process PUSH_TAXA_FILES {
    tag "${taxa}"
    label 'process_low'

    input:
    tuple val(taxa), path(db), val(db_type)

    output:
    tuple path(db), val(db_type), emit: taxa_files

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """
}
