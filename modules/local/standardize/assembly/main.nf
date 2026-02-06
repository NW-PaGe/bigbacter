process STD_ASSEMBLY {
    tag "${meta.id}"
    label "process_low"
    stageInMode "copy"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.fa.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[ "${assembly}" == *.gz ]]; then
        # Already gzipped: just rename/move
        cp "${assembly}" "${meta.id}.fa.gz"
    else
        pigz -p "${task.cpus}" "${assembly}"
        cp "${assembly}.gz" "${meta.id}.fa.gz"
    fi
    """
}