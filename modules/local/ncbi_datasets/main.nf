process NCBI_DATASETS {
    label 'process_low'

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.fna"), emit: assembly
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    """
    # download assembly, unzip, and move for publishing
    datasets download genome accession ${accession}
    unzip ncbi_dataset.zip
    mv ncbi_dataset/data/*/*.fna ./

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets: \$(datasets --version | sed 's/.*: //g')
    END_VERSIONS
    """
}