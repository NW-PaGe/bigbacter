process PAN_DISTS {
    label "process_low"

    input:
    tuple val(meta), path(panaroo, stageAs: 'panaroo/*')

    output:
    tuple val(meta), path("pan-dists/*"),                      emit: results
    tuple val(meta), path("pan-dists/dist_events_all.tsv"),    emit: dist_events_all
    tuple val(meta), path("pan-dists/dist_events_multi.tsv"),  emit: dist_events_multi
    tuple val(meta), path("pan-dists/dist_genes.tsv"),         emit: dist_genes
    tuple val(meta), path("pan-dists/pairwise.tsv"),           emit: pairwise
    tuple val(meta), path("pan-dists/events.tsv"),             emit: events
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    tool = 'pan-dists.py'
    """
    # Accessory gain/loss events from the Panaroo pangenome, and the pairwise
    # distances derived from them. Reads gene_presence_absence_roary.csv for
    # cluster membership and gene_data.csv for gene order along each contig.
    ${tool} \\
        panaroo \\
        -o pan-dists \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """

    stub:
    """
    mkdir pan-dists
    printf 'event_id\\tevent_type\\tn_clusters\\n' > pan-dists/events.tsv
    printf 'isolate_a\\tisolate_b\\tevents_all\\tevents_multi\\tevents_single\\n' > pan-dists/pairwise.tsv
    for f in dist_events_all dist_events_multi dist_genes; do
        printf '\\n' > pan-dists/\$f.tsv
    done
    touch pan-dists/detail.tsv pan-dists/summary.txt pan-dists/report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: stub
    END_VERSIONS
    """
}
