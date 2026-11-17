process BAKTA {
    label "process_medium"

    input:
    tuple val(meta), path(assembly), path(db)

    output:
    tuple val(meta), path("${prefix}.gff.gz"), emit: gff
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = "${meta.id}"
    """
    bakta \\
        --db ${db} \\
        --output bakta \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        --force \\
        ${args} \\
        ${assembly}

    # only the GFF3 (annotations + sequences) is kept
    gzip -c bakta/${prefix}.gff3 > ${prefix}.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: "\$(bakta --version 2>/dev/null | sed -n 's/^bakta //p' | tail -n 1)"
    END_VERSIONS
    """

    stub:
    prefix = "${meta.id}"
    """
    echo "##gff-version 3" | gzip -c > ${prefix}.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: stub
    END_VERSIONS
    """
}
