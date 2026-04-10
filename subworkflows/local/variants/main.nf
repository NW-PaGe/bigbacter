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
            [meta, ref_path, ref_meta, []]
        }

    ch_cluster_data_old = ch_cluster_status
        .old_cluster
        .map { meta ->
            def cluster_path = file(params.db) / meta.taxa / 'clusters' / meta.cluster
            def ref_path     = cluster_path / 'asm' / 'ref.fa.gz'
            def ref_meta     = cluster_path / 'aux' / 'ref.json'
            def snp_dir      = cluster_path / 'var'

            if (!ref_path.exists()) {
                error "Reference assembly not found for cluster ${meta.cluster}: ${ref_path}"
            }

            def ref_meta_file = ref_meta.exists() ? ref_meta : []
            def snp_files     = snp_dir.exists()
                ? snp_dir.listFiles().findAll { it.name.endsWith('.tar.gz') }
                : []

            [meta, ref_path, ref_meta_file, snp_files]
        }

    ch_cluster_data = ch_cluster_data_new.concat(ch_cluster_data_old)

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
            .map { [Utils.cluster_meta(it.meta), it.meta, it.reads, it.assembly] }
            .combine(ch_cluster_data, by: 0)
            .map { cluster_meta, meta, reads, assembly, ref_path, ref_meta, snp_files ->
                [meta.id, meta, reads, assembly, ref_path]
            }
            .join(ch_ds_rate, by: 0, remainder: true)
            .map { sid, meta, reads, assembly, ref_path, rate ->
                [meta, reads[0] ? reads : [], assembly, ref_path, rate]
            }
    )
    ch_versions = ch_versions.mix(SNIPPY_SINGLE.out.versions.first())

    /*
    =============================================================================================================================
        CLUSTER-LEVEL CORE SNP ALIGNMENT
    =============================================================================================================================
    */

    ch_indv_snp_files = SNIPPY_SINGLE
        .out
        .results
        .map { meta, new_snps -> [Utils.cluster_meta(meta), new_snps] }
        .groupTuple()
        .join(ch_cluster_data)
        .map { meta, new_snps, ref_path, ref_meta, old_snps ->
            def new_names     = new_snps.collect { it.name }
            def old_snps_keep = old_snps.flatten().findAll { !new_names.contains(it.name) }
            def all_snps      = (new_snps + old_snps_keep).flatten()

            [meta + [n: all_snps.size()], ref_path, all_snps]
        }

    SNIPPY_CORE(ch_indv_snp_files)
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions.first())

    ch_aln = SNIPPY_CORE
        .out
        .full_aln
        .map { meta, aln ->
            [meta + [recomb_masked: false], aln, []]
        }

    /*
    =============================================================================================================================
        RECOMBINATION MASKING (OPTIONAL)
    =============================================================================================================================
    */

    if (params.mask_recomb) {

        ch_aln_clean_filt = SNIPPY_CORE
            .out
            .clean_aln
            .filter { meta, aln -> meta.n > params.min_tree }

        GUBBINS(ch_aln_clean_filt)
        ch_versions = ch_versions.mix(GUBBINS.out.versions.first())

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
    snp_files   = SNIPPY_SINGLE.out.results
    snp_summary = POLYCORE.out.csv
    snp_dist    = POLYCORE.out.dist_wide
        .join(POLYCORE.out.dist_long)
        .map { meta, dist_wide, dist_long -> [meta, [dist_wide, dist_long]] }
    snp_tree    = IQTREE.out.result
    versions    = ch_versions
}
