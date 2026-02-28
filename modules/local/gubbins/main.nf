process GUBBINS {
    tag "${taxa}-${cluster}"
    label 'process_medium'
    errorStrategy 'ignore'

    conda "bioconda::gubbins=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins%3A3.3.1--py39pl5321h3d4b85c_0' :
        'quay.io/biocontainers/gubbins:3.3.1--py310pl5321h83093d7_0' }"

    containerOptions "${ workflow.containerEngine != 'apptainer' && workflow.containerEngine != 'singularity' ? '--shm-size 1000000000' : '' }"

    input:
    tuple val(taxa), val(cluster), path(aln), path(const_sites), val(count)
    val timestamp

    output:
    tuple val(taxa), val(cluster), path("*.fasta"), path(const_sites), val(count), emit: aln
    tuple val(taxa), val(cluster), path("*.gff"),                                  emit: gff
    tuple val(taxa), val(cluster), path("*.vcf"),                                  emit: vcf
    tuple val(taxa), val(cluster), path("*.gubbins.stats"),                        emit: stats
    tuple val(taxa), val(cluster), path("*.recombination_predictions.embl"),       emit: embl_predicted
    tuple val(taxa), val(cluster), path("*.branch_base_reconstruction.embl"),      emit: embl_branch
    path "versions.yml",                                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = "${timestamp}-${taxa}-${cluster}"
    """
    # set resource limit
    ulimit -m ${task.memory.toBytes()}

    # Run Gubbins with --no-cleanup to prevent temporary file deletion
    run_gubbins.py \\
        --threads $task.cpus \\
        --prefix ${prefix} \\
        --no-cleanup \\
        ${aln}

    # rename stats for easy summary
    mv *.per_branch_statistics.csv ${prefix}.gubbins.stats
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1)
    END_VERSIONS
    """
}
