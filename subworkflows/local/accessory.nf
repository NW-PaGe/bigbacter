//
// Perform accessory genome analysis within each cluster
//

include { FORMAT_DIST } from '../../modules/local/format-dist'
include { DIST_MAT    } from '../../modules/local/dist-mat'

workflow ACCESSORY {
    take:
    core_acc_dist // channel: [ val(taxa), val(dist) ]
    core_tree     // channel: [ val(taxa), val(cluster), path(tree), val(type), val(method), val(source), val(stats)]
    manifest_file // channel: [ val(new_samples) ]
    timestamp     // channel: val(timestamp)

    main:
    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        SUBSET & FORMAT ACCESSORY GENOME DISTANCES
        - Bash is used for subsetting because it is much faster than Python or R
    =============================================================================================================================
    */ 
    // Select tree generated using Snippy alignment (ML or NJ)
    core_tree
        .filter{ taxa, cluster, source, tree -> source == "snippy" }
        .map{ taxa, cluster, source, tree -> [ taxa, cluster, tree ] }
        .set{ core_tree }
    // Combine tree with the full PopPUNK distance file
    core_acc_dist
        .combine(core_tree, by: 0)
        .map{ taxa, dist, cluster, tree -> [ taxa, cluster, dist, tree ] }
        .set { full_dist }
    // MODULE: Subset and format dist file using Snippy tree and Bash
    FORMAT_DIST (
        full_dist
    )
    
    /* 
    =============================================================================================================================
        CREATE ACCESSORY DISTANCE MATRIX FIGURE
    =============================================================================================================================
    */ 
    // MODULE: Distance matrix figure
    DIST_MAT (
        FORMAT_DIST.out.dist.map{ taxa, cluster, dist, tree -> [ taxa, cluster, dist.baseName.endsWith('.csv') ? 'bb-cluster' : "poppunk", dist, tree ] }.combine(manifest_file),
        "long",
        "Accessory",
        100,
        "TRUE",
        timestamp
    )

    emit:
    dist      = DIST_MAT.out.dist_wide // channel: [ val(taxa), val(cluster), val(input_source), path(dist) ] 
    versions  = ch_versions            // channel: [ versions.yml ]
}
