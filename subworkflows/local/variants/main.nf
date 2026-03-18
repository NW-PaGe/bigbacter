//
// Core genome SNP analysis and phylogenetic reconstruction
//

include { SELECT_REF    } from '../../../modules/local/select_ref/main'
include { DS_RATE       } from '../../../modules/local/ds_rate/main'
include { SNIPPY_SINGLE } from '../../../modules/local/snippy/single/main'
include { SNIPPY_CORE   } from '../../../modules/local/snippy/core/main'
include { POLYCORE      } from '../../../modules/local/polycore/main'
include { GUBBINS       } from '../../../modules/local/gubbins/main'
include { IQTREE        } from '../../../modules/local/iqtree/main'

/*
=============================================================================================================================
    SUBWORKFLOW
=============================================================================================================================
*/

workflow VARIANTS {

    take:
    ch_manifest   // channel: [ val(meta), path(assembly), path(reads) ]
    ch_read_stats // channel: [ val(meta), path(fastp_json) ]

    main:
    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        IDENTIFY NEW VS EXISTING CLUSTERS
    =============================================================================================================================
    */

    ch_cluster_status = ch_manifest
        .map { Utils.cluster_meta(it.meta) }
        .unique()
        .branch { meta ->
            def cluster_path   = file(params.db) / meta.taxa / 'clusters' / meta.cluster
            def cluster_exists = cluster_path.exists()

            new_cluster: !cluster_exists
            old_cluster: cluster_exists
        }
        
    /*
    =============================================================================================================================
        PREPARE CLUSTER DATA (REFERENCE + PREVIOUS SNPS)
    =============================================================================================================================
    */

    // Select reference genome for new clusters
    SELECT_REF(
        ch_manifest
            .map { [Utils.cluster_meta(it.meta), it.assembly] }
            .groupTuple(by: 0)
            .join(ch_cluster_status.new_cluster)
    )
    ch_versions = ch_versions.mix(SELECT_REF.out.versions.first())

    ch_cluster_data_new = SELECT_REF
        .out
        .result
        .map { meta, ref_path, ref_meta -> 
            [meta, ref_path, ref_meta, []] // empty list for snp_files
        }

    // Gather reference and SNP files for existing clusters
    ch_cluster_data_old = ch_cluster_status
        .old_cluster
        .map { meta ->
            def cluster_path = file(params.db) / meta.taxa / 'clusters' / meta.cluster
            def ref_path     = cluster_path / 'asm' / 'ref.fa.gz'
            def ref_meta     = cluster_path / 'aux' / 'ref.json'
            def snp_dir      = cluster_path / 'var'

            // Validate reference exists
            if (!ref_path.exists()) {
                error "Reference assembly not found for cluster ${meta.cluster}: ${ref_path}"
            }

            // Load reference metadata if available
            def ref_meta_file = ref_meta.exists() ? ref_meta : []

            // Gather existing SNP tarballs
            def snp_files = snp_dir.exists() 
                ? snp_dir.listFiles().findAll { it.name.endsWith('.tar.gz') }
                : []

            [meta, ref_path, ref_meta_file, snp_files]
        }

    ch_cluster_data = ch_cluster_data_new.concat(ch_cluster_data_old)

    /*
    =============================================================================================================================
        PRESERVE TIMESTAMP FOR ALL CLUSTER-LEVEL PROCESSES
    =============================================================================================================================
    */

    // Create a reusable channel mapping cluster -> timestamp
    // This will be joined with ALL cluster-level metadata to ensure publishDir paths are correct
    ch_cluster_timestamps = ch_manifest
        .map { tuple ->
            def cluster_meta = Utils.cluster_meta(tuple.meta)
            def timestamp = tuple.meta.timestamp
            [[taxa: cluster_meta.taxa, cluster: cluster_meta.cluster], timestamp]
        }
        .unique()

    /*
    =============================================================================================================================
        CALCULATE READ DOWNSAMPLING RATES
    =============================================================================================================================
    */

    DS_RATE(
        ch_read_stats
            .map { meta, json -> [meta.id, json] }
            .join(
                ch_manifest.map { [it.meta.id, Utils.cluster_meta(it.meta)] }
            )
            .map { sid, json, cluster_meta -> [cluster_meta, json] }
            .groupTuple()
            .map { cluster_meta, jsons ->
                // Create join key from cluster_meta
                [[taxa: cluster_meta.taxa, cluster: cluster_meta.cluster], cluster_meta, jsons]
            }
            .join(ch_cluster_timestamps)
            .map { key, cluster_meta, jsons, timestamp ->
                // Add timestamp back to cluster_meta
                [cluster_meta + [timestamp: timestamp], jsons]
            }
            .join(
                ch_cluster_data.map { meta, ref_path, ref_meta, snp_files -> 
                    [Utils.cluster_meta(meta), ref_path] 
                }
            )
    )
    ch_versions = ch_versions.mix(DS_RATE.out.versions.first())

    ch_ds_rate = DS_RATE
        .out
        .result
        .map { meta, csv -> csv }
        .splitCsv(header: true)
        .map { [it.filename, it.rate] }
        .join(ch_read_stats.map { meta, json -> [json.name, meta] })
        .map { filename, rate, meta -> [meta.id, rate] }
        
    /*
    =============================================================================================================================
        PER-SAMPLE SNP CALLING
    =============================================================================================================================
    */

    SNIPPY_SINGLE(
        ch_manifest
            .map { [Utils.cluster_meta(it.meta), it.meta, it.reads] }
            .combine(ch_cluster_data, by: 0)
            .map { cluster_meta, meta, reads, ref_path, ref_meta, snp_files -> 
                [meta.id, meta, reads, ref_path] 
            }
            .join(ch_ds_rate, by: 0)
            .map { sid, meta, reads, ref_path, rate -> 
                [meta, reads, ref_path, rate] 
            }
    )
    ch_versions = ch_versions.mix(SNIPPY_SINGLE.out.versions.first())

    /*
    =============================================================================================================================
        CLUSTER-LEVEL CORE SNP ALIGNMENT
    =============================================================================================================================
    */

    // Combine new and existing SNP calls per cluster
    ch_indv_snp_files = SNIPPY_SINGLE
        .out
        .results
        .map { meta, new_snps -> [Utils.cluster_meta(meta), new_snps] }
        .groupTuple()
        .join(ch_cluster_data)                    
        .map { meta, new_snps, ref_path, ref_meta, old_snps ->
            // Exclude old SNP files that have been regenerated
            def new_names     = new_snps.collect { it.name }
            def old_snps_keep = old_snps.flatten().findAll { !new_names.contains(it.name) }
            def all_snps      = (new_snps + old_snps_keep).flatten()

            // Create join key for timestamp enrichment
            [[taxa: meta.taxa, cluster: meta.cluster], meta + [n: all_snps.size()], ref_path, all_snps]
        }
        .join(ch_cluster_timestamps)
        .map { key, meta, ref_path, all_snps, timestamp ->
            // Add timestamp to metadata before SNIPPY_CORE
            [meta + [timestamp: timestamp], ref_path, all_snps]
        }

    SNIPPY_CORE(ch_indv_snp_files)
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions.first())

    // Prepare alignment channel with recombination mask placeholder
    // Note: meta already has timestamp from SNIPPY_CORE enrichment above
    ch_aln = SNIPPY_CORE
        .out
        .full_aln
        .map { meta, aln -> 
            [meta + [recomb_masked: false], aln, []] // empty mask file
        }

    /*
    =============================================================================================================================
        RECOMBINATION MASKING (OPTIONAL)
    =============================================================================================================================
    */

    if (params.mask_recomb) {
        
        // Filter clusters with sufficient samples for recombination detection
        ch_aln_clean_filt = SNIPPY_CORE
            .out
            .clean_aln
            .filter { meta, aln -> meta.n > params.min_tree }

        // Detect recombinant regions with Gubbins
        GUBBINS(ch_aln_clean_filt)
        ch_versions = ch_versions.mix(GUBBINS.out.versions.first())

        // Add recombination-masked alignments to channel
        ch_aln = ch_aln_clean_filt
            .join(GUBBINS.out.bed)
            .map { meta, aln, mask -> 
                [meta + [recomb_masked: true], aln, mask] 
            }
            .concat(ch_aln)
    }

    /*
    =============================================================================================================================
        CORE GENOME STATISTICS
    =============================================================================================================================
    */

    POLYCORE(ch_aln)
    ch_versions = ch_versions.mix(POLYCORE.out.versions.first())

    /*
    =============================================================================================================================
        PHYLOGENETIC RECONSTRUCTION
    =============================================================================================================================
    */

    // Build ML trees for clusters with sufficient samples
    IQTREE(
        POLYCORE
            .out
            .aln
            .join(POLYCORE.out.fconst)
            .filter { meta, aln, fconst -> meta.n > params.min_tree }
    )
    ch_versions = ch_versions.mix(IQTREE.out.versions.first())

    emit:
    ref_data    = ch_cluster_data.map { meta, ref_path, ref_meta, snp_files -> 
        [meta, ref_path, ref_meta] 
    }
    snp_files   = SNIPPY_SINGLE.out.results   // channel: [ val(meta), path(snp_tar) ]
    snp_summary = POLYCORE.out.csv            // channel: [ val(meta), path(summary_csv) ]
    snp_dist    = POLYCORE.out.dist_wide      // channel: [ val(meta), path(distance_matrix) ]
    snp_tree    = IQTREE.out.result           // channel: [ val(meta), path(tree) ]
    versions    = ch_versions                 // channel: [ path(versions.yml) ]
}
