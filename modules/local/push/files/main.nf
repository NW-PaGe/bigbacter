process FILES {
    input:
    tuple val(meta), path(files)

    output:
    path "*", includeInputs: true

    script:
    """
    for f in *; do
        [ -e "\$f" ] || continue
        pigz -p ${task.cpus} "\$f" || true
    done
    """
}