
process RESOLVE_MERGED {
    tag "${sample}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(sample), path(dist)

    output:
    tuple val(taxa), val(sample), env(NEIGHBOR), emit: neighbor

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    zcat ${dist} \\
        | awk -v OFS='\t' -v sample=${sample} '! (\$1 == sample && \$2 == sample) && (\$1 == sample | \$2 == sample) {print \$1, \$2, \$5}' \\
        | sort -rgk 2 \\
        | sed -n 1p \\
        | cut -f 1,2 \\
        | tr '\t' '\n' \\
        | grep -v ${sample} \\
        | sort \\
        | uniq \\
        | read NEIGHBOR
    """
}