process FASTERQDUMP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), val(sra)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fasterq-dump \\
        $args \\
        --threads $task.cpus \\
        ${sra}

    pigz \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}