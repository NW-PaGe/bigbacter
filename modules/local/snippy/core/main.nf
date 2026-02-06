process SNIPPY_CORE {
    tag "${meta.taxa}-${meta.cluster}"
    label "process_high"
    stageInMode "copy"

    input:
    tuple val(meta), path(ref, stageAs: 'ref.fa.gz'), path(snp_files, stageAs: 'snps_tar/*')

    output:
    tuple val(meta), path("${prefix}.aln"),       emit: aln,       optional: true
    tuple val(meta), path("${prefix}.full.aln"),  emit: full_aln,  optional: true
    tuple val(meta), path("${prefix}.clean.aln"), emit: clean_aln, optional: true
    tuple val(meta), path("${prefix}.txt"),       emit: txt
    tuple val(meta), path("${prefix}.tab"),       emit: tab
    tuple val(meta), path("${prefix}.vcf"),       emit: vcf
    path 'versions.yml',                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = "${meta.timestamp}-${meta.taxa}-${meta.cluster}"
    """
    gzip -d ${ref}

    # extract file contents
    mkdir snps_extracted
    for FILE in ${snp_files}; do
        echo "Extracting \${FILE}"
        tar -xzf \${FILE} -C snps_extracted/
    done

    # run Snippy Core - fails when no SNPs are found but we let it pass
    snippy-core \\
        --prefix "${prefix}" \\
        --ref ref.fa \\
        ${args} \\
        snps_extracted/*/ \\
        || true

    # create clean version for Gubbins
    snippy-clean_full_aln ${prefix}.full.aln > ${prefix}.clean.aln
    
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snippy: \$(snippy --version | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
