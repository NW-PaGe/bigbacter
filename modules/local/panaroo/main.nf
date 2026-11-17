process PANAROO {
    label "process_high"

    input:
    tuple val(meta), path(gffs, stageAs: 'gff_gz/*')

    output:
    tuple val(meta), path("panaroo/*"),                           emit: results
    tuple val(meta), path("panaroo/gene_presence_absence.Rtab"),  emit: rtab
    tuple val(meta), path("panaroo/gene_presence_absence.csv"),   emit: csv
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Panaroo does not read gzipped GFFs. Samples are named after the files (<sample>.gff).
    mkdir gff
    for f in gff_gz/*.gff.gz; do
        gzip -dc "\$f" > "gff/\$(basename "\${f%.gz}")"
    done

    panaroo \\
        -i gff/*.gff \\
        -o panaroo \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: "\$(panaroo --version 2>/dev/null | sed -n 's/^panaroo //p' | tail -n 1)"
    END_VERSIONS
    """

    stub:
    """
    mkdir panaroo
    ls gff_gz | sed 's/.gff.gz\$//' > panaroo/samples.txt
    touch panaroo/gene_presence_absence.Rtab panaroo/gene_presence_absence.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: stub
    END_VERSIONS
    """
}
