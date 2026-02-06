process SNIPPY_SINGLE {
    tag "${meta.id}"
    label "process_high"
    stageInMode "copy"

    input:
    tuple val(meta), path(reads), path(ref, stageAs: 'ref.fa.gz'), val(downsample_rate)

    output:
    tuple val(meta), path('*.tar.gz'), emit: results

    path 'versions.yml',               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    """
    # extract reference
    gzip -d ${ref}

    # run Snippy
    snippy \\
        --reference ref.fa \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        --outdir ${meta.id} \\
        --cpus ${task.cpus} \\
        --subsample ${downsample_rate} \\
        ${args}

    # compress output
    tar -czvf ${meta.id}.tar.gz ${meta.id}/

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snippy: \$(snippy --version | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
