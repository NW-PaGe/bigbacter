process GUBBINS {
    label 'process_medium'
    errorStrategy 'ignore'

    containerOptions "${ workflow.containerEngine != 'apptainer' && workflow.containerEngine != 'singularity' ? '--shm-size 1000000000' : '' }"

    input:
    tuple val(meta), path(aln)

    output:
    tuple val(meta), path("*.fasta"),                                emit: aln
    tuple val(meta), path("*.gff"),                                  emit: gff
    tuple val(meta), path("*.bed"),                                  emit: bed
    tuple val(meta), path("*.vcf"),                                  emit: vcf
    tuple val(meta), path("*.per_branch_statistics.csv"),            emit: stats
    tuple val(meta), path("*.recombination_predictions.embl"),       emit: embl_predicted
    tuple val(meta), path("*.branch_base_reconstruction.embl"),      emit: embl_branch
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # set resource limit
    ulimit -m ${task.memory.toBytes()}

    # Run Gubbins
    run_gubbins.py \\
        --threads $task.cpus \\
        ${aln}

    # convert gff to bed
    cat *.gff | grep -v '#' | awk '{print "Reference\\t"\$4"\\t"\$5}' > gubbins.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2> /dev/null)
    END_VERSIONS
    """
}