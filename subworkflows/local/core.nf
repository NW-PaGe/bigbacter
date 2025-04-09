//
// Perform core genome analysis within each cluster
//

include { SNIPPY_SINGLE } from '../../modules/local/snippy'
include { SNIPPY_CORE   } from '../../modules/local/snippy'
include { GUBBINS       } from '../../modules/local/gubbins'
include { SNP_DISTS     } from '../../modules/local/snp-dists'
include { IQTREE        } from '../../modules/local/iqtree'
include { RAPIDNJ       } from '../../modules/local/rapidnj'
include { TREE_FIGURE   } from '../../modules/local/tree-figures'
include { DIST_MAT      } from '../../modules/local/dist-mat'

/*
=============================================================================================================================
    SUBWORKFLOW FUNCTIONS
=============================================================================================================================
*/

// Function for counting the number of samples in an alignment file
def count_alignments ( aln_file ) {
    count = 0
    total = aln_file.eachLine{ line -> count+= line.count('>')}
    return total
}

/*
=============================================================================================================================
    SUBWORKFLOW
=============================================================================================================================
*/
workflow CORE {
    take:
    ch_man      // channel: [ val(sample), val(taxa), path(assembly), path(fastq_1), path(fastq_2), val(cluster) ]
    ch_man_file // channel: [ val(samplesheet.valid.csv) ]
    timestamp     // channel: val(timestamp)

    main:
    ch_versions = Channel.empty()
    /* 
    =============================================================================================================================
        DETERMINE REFERENCE GENOME
    =============================================================================================================================
    */
    // Group by taxa and cluster
    ch_man
        .map{ sample, taxa, assembly, fastq_1, fastq_2, cluster -> [ taxa, cluster, assembly ]}
        .groupTuple(by: [0,1])
        .map{ taxa, cluster, assembly -> [ taxa, cluster, assembly, file(params.db).resolve(taxa).resolve("clusters").resolve(cluster).exists() ? true : false ] }
        .set { clust_grps }
    // New clusters - select first assembly
    clust_grps
        .filter{ taxa, cluster, assembly, status -> ! status }
        .map{ taxa, cluster, assembly, status -> [ taxa, cluster, assembly.first(), status ] }
        .set{ clust_grp_new }
    // Old clusters - resolve path to reference in BigBacter database
    clust_grps
        .filter{ taxa, cluster, assembly, status -> status }
        .map{ taxa, cluster, assembly, status -> [ taxa, cluster, file(params.db) / taxa / "clusters" / cluster / "ref/ref.fa.gz", status ] }
        .set{ clust_grps_old }
    // Combine new and old clusters into single channel and add back to reference
    ch_man
        .map { sample, taxa, assembly, fastq_1, fastq_2, cluster -> [ taxa, cluster, sample, assembly, fastq_1, fastq_2 ] }
        .combine( clust_grp_new.concat( clust_grps_old ), by: [0,1] )
        .map { taxa, cluster, sample, assembly, fastq_1, fastq_2, ref, status -> [sample, taxa, assembly, fastq_1, fastq_2, cluster, status, ref] }
        .set{ ch_man }

    /* 
    =============================================================================================================================
        INDIVIDUAL SNP CALLING
    =============================================================================================================================
    */
    // Run Snippy on each isolate
    SNIPPY_SINGLE(
        ch_man,
        timestamp
    )
    ch_versions = ch_versions.mix(SNIPPY_SINGLE.out.versions)

    // Group samples by taxa and cluster
    ch_man
        .map{ sample, taxa, assembly, fastq_1, fastq_2, cluster, status, ref -> [ sample, taxa, cluster, status, ref ] }
        .join(SNIPPY_SINGLE.out.results, by: [0,1])
        .groupTuple(by: [1,2])
        .map {sample, taxa, cluster, status, ref, new_snps -> [ taxa, cluster, status.get(0), ref.get(0), new_snps ]}
        .set {clust_grps}
    // New clusters - add empty tuple in place of historic SNP files
    clust_grps
        .filter{ taxa, cluster, status, ref, new_snps -> ! status }
        .map{ taxa, cluster, status, ref, new_snps -> [ taxa, cluster, ref, new_snps, [] ] }
        .set{ clust_grp_new }
    // Old clusters - resolve path to snippy directory in the BigBacter database
    clust_grps
        .filter{ taxa, cluster, status, ref, new_snps -> status }
        .map{ taxa, cluster, status, ref, new_snps -> [ taxa, cluster, ref, new_snps, file(file(params.db) / taxa / "clusters" / cluster / "snippy", type: 'dir') ] }
        .set{ clust_grp_old }
    // Combine new and old cluster channels
    clust_grp_new
        .concat(clust_grp_old)
        .set { snp_files }
    
    /* 
    =============================================================================================================================
        CORE SNP CALLING
    =============================================================================================================================
    */
    // Run Snippy Core
    SNIPPY_CORE(
        snp_files,
        timestamp
    )
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)

    /* 
    =============================================================================================================================
        GATHER CORE ALIGNMENT METRICS
    =============================================================================================================================
    */
    // count the number of samples in each alignment & determine if any SNPs could be identified
    SNIPPY_CORE
        .out
        .aln
        .map{ taxa, cluster, core_aln, clean_aln, const_sites -> [ taxa, cluster, core_aln, clean_aln, const_sites, count_alignments(core_aln), file(core_aln).isEmpty() ] }
        .filter{ taxa, cluster, core_aln, clean_aln, const_sites, count, snp_status -> ! snp_status && count > params.min_tree } // filter alignments with no SNPs and those below the min number of samples
        .map{ taxa, cluster, core_aln, clean_aln, const_sites, count, snp_status -> [ taxa, cluster, core_aln, clean_aln, const_sites, count ] }
        .set{ aln_w_metrics }

    

       /*
    =============================================================================================================================
        IDENTIFY RECOMBINANT REGIONS
        - only ran for clusters smaller than '--max_ml'
        - alignments without SNPs are excluded
    =============================================================================================================================
    */
    // MODULE: Run Gubbins
    GUBBINS(
        aln_w_metrics.filter{ taxa, cluster, core_aln, clean_aln, const_sites, count -> count <= params.max_ml }.map{ taxa, cluster, core_aln, clean_aln, const_sites, count -> [ taxa, cluster, clean_aln, const_sites, count ] },
        timestamp
    )
    ch_versions = ch_versions.mix(GUBBINS.out.versions)

    // Add Gubbins alignments to channel
    aln_w_metrics
        .map{ taxa, cluster, core_aln, clean_aln, const_sites, count -> [ taxa, cluster, core_aln, const_sites, count, "snippy" ]} // no longer need the full alignments
        .concat(GUBBINS.out.aln.map{ taxa, cluster, core_aln, const_sites, count -> [ taxa, cluster, core_aln, const_sites, count, "gubbins" ] })
        .set{ aln_w_metrics }
    /*
    =============================================================================================================================
        GENERATE NEIGHBOR JOINING TREE
        - only performed for clusters with SNPs, samples > 'max_ml',and samples > 'min_tree'
        - this option should only be used if you do not have the computational power needed for a ML tree
    =============================================================================================================================
    */
    // Filter clusters exceeding max_ml and above min_tree
    aln_w_metrics
        .filter{ taxa, cluster, aln, const_sites, count, source -> count > params.max_ml }
        .map{ taxa, cluster, aln, const_sites, count, source -> [ taxa, cluster, aln, source ]}
        .set{ nj_samples }
    // MODULE: Run Rapidnj 
    RAPIDNJ(
        nj_samples,
        timestamp
    )
    ch_versions = ch_versions.mix(RAPIDNJ.out.versions)
    /* 
    =============================================================================================================================
        GENERATE MAXIMUM LIKELIHOOD TREE
        - only performed for clusters with SNPs, samples <= 'max_ml', and samples > 'min_tree'
        - this is the preferred method
        - this allows for bootstrapping of Snippy tree - which is not performed on initial iterations by Gubbins
    =============================================================================================================================
    */
    // Filter clusters at or below max_ml and above min_tree
    aln_w_metrics
        .filter{ taxa, cluster, aln, const_sites, count, source -> count <= params.max_ml }
        .map{ taxa, cluster, aln, const_sites, count, source -> [ taxa, cluster, aln, const_sites, count, source ]}
        .set{ ml_samples }
    // MODULE: Run IQTREE2
    IQTREE(
        ml_samples,
        timestamp
    )
    ch_versions = ch_versions.mix(IQTREE.out.versions)

    /* 
    =============================================================================================================================
        CALCULATE PAIRWISE SNP DIFFERENCES
        - Performed on full Snippy alignments and core Gubbins alignments
    =============================================================================================================================
    */
    // Combine alignment files
    SNIPPY_CORE
        .out
        .aln
        .map{ taxa, cluster, core_aln, clean_aln, const_sites -> [ taxa, cluster, clean_aln, "snippy" ] }
        .concat( GUBBINS.out.aln.map{ taxa, cluster, core_aln, const_sites, count -> [ taxa, cluster, core_aln, "gubbins" ] })
        .set{ all_alns }

    // Create SNP distance matrix
    SNP_DISTS(
        all_alns,
        timestamp
    )
    ch_versions = ch_versions.mix(SNP_DISTS.out.versions)

    /*
    =============================================================================================================================
        CREATE TREE FIGURES
    =============================================================================================================================
    */
    // Trees Channel
    RAPIDNJ
        .out
        .result
        .map{ taxa, cluster, tree, source -> [taxa, cluster, tree, source, "Core SNPs", "Neighbor Joining", [] ] }
        .set{ nj_tree }
    IQTREE
        .out
        .result
        .map{ taxa, cluster, tree, source -> [taxa, cluster, tree, source, "Core SNPs", "Maximum Likelihood" ] }
        .combine(SNIPPY_CORE.out.stats, by: [0,1])
        .set{ ml_tree }
    ml_tree
        .concat(nj_tree)
        .set{ all_trees }
    
    // Tree figures
    TREE_FIGURE (
        all_trees.combine(ch_man_file),
        timestamp
    )

    /* 
    =============================================================================================================================
        CREATE SNP DISTANCE MATRIX FIGURES
        - boolean indicates if recombination is masked
    =============================================================================================================================
    */
    // Combine distance matrix with tree
    SNP_DISTS
        .out
        .result
        .map{ taxa, cluster, dist, source -> [ taxa, cluster, source, dist ] }
        .join(all_trees.map{ taxa, cluster, tree, source, type, method, stats -> [ taxa, cluster, source, tree ]}, by: [0,1,2], remainder: true)
        .map{ taxa, cluster, source, dist, tree -> [ taxa, cluster, source, dist, tree == null ? [] : tree ] }       
        .set{ all_dists }

    // Create distance matrix figure
    DIST_MAT (
        all_dists.combine(ch_man_file),
        "wide",
        "Core SNPs",
        100,
        "FALSE",
        timestamp
    )

    /* 
    =============================================================================================================================
        COMBINE CHANNEL FOR EMITTING
    =============================================================================================================================
    */
    // Combine stats files from core SNP analysis
    SNIPPY_CORE
        .out
        .stats
        .concat(GUBBINS.out.stats)
        .groupTuple(by: [0,1])
        .set{ all_stats }

    emit:
    snp_files = snp_files              // channel: [ val(taxa), val(cluster), path(ref), path(new_snippy), path(old_snippy) ]
    tree      = TREE_FIGURE.out.tree   // channel: [ val(taxa), val(cluster), val(source), path(tree) ]
    dist      = DIST_MAT.out.dist_wide // channel: [ val(taxa), val(cluster), val(source), path(dist) ]
    meta      = TREE_FIGURE.out.meta      // channel: [ val(taxa), val(cluster), val(source), path(meta) ]
    stats     = all_stats              // channel: [ val(taxa), val(cluster), path(stats) ]
    versions  = ch_versions            // channel: [ versions.yml ]
}
