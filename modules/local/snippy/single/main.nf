process SNIPPY_SINGLE {
    tag "${meta.id}"
    label "process_high"
    stageInMode "copy"

    input:
    tuple val(meta), path(reads), path(assembly), path(ref, stageAs: 'ref.fa.gz'), val(downsample_rate)

    output:
    tuple val(meta), path('*.tar.gz'), emit: results

    path 'versions.yml',               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    asm_gz_arg    = reads ? "" : "gzip -d ${assembly}"
    reads_arg     = reads ? "--R1 ${reads[0]} --R2 ${reads[1]}" : "--ctgs ${assembly.baseName}"
    subsample_arg = reads ? "--subsample ${downsample_rate}" : ""
    cleanup_arg   = params.keep_bam ? "" : "--cleanup"
    
    """
    # extract reference
    gzip -d ${ref}

    # extract assembly (if used)
    ${asm_gz_arg}

    # run Snippy
    snippy \\
        --reference ref.fa \\
        ${reads_arg} \\
        ${subsample_arg} \\
        --outdir ${meta.id} \\
        --cpus ${task.cpus} \\
        ${cleanup_arg} \\
        ${args}

    # compress output
    tar -czvf ${meta.id}.tar.gz ${meta.id}/

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snippy: \$(snippy --version | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
