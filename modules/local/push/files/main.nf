process FILES {
    input:
    tuple val(meta), path(files)

    output:
    path "*.fa.gz",   includeInputs: true, optional: true
    path "*.tar.gz",  includeInputs: true, optional: true
    path "*.json.gz", includeInputs: true, optional: true
    path "*.sig.gz",  includeInputs: true, optional: true

    script:
    """
    for f in *; do
        [ -e "\$f" ] || continue
        pigz -p ${task.cpus} "\$f" || true
    done
    """
}