//
// Taxonomic classification and cluster assignment workflow
//

include { GAMBIT_QUERY } from '../../../modules/local/gambit/main'
include { FLOC         } from '../../../modules/local/floc/main'

workflow CLASSIFY {
    take:
    ch_manifest // channel: [ val(meta), path(assembly), ... ]
    
    main:
    ch_versions = Channel.empty()

    /* 
    =============================================================================================================================
        SPECIES PREDICTION (GAMBIT)
    =============================================================================================================================
    */

    // Predict taxa for samples without pre-assigned taxonomy
    GAMBIT_QUERY(
        ch_manifest
            .filter { !it.meta.taxa }
            .map { [it.meta, it.assembly] },
        params.gambit_db
    )
    ch_versions = ch_versions.mix(GAMBIT_QUERY.out.versions.first())
    
    // Parse GAMBIT results and enrich metadata with predicted taxa
    ch_manifest_with_taxa = GAMBIT_QUERY
        .out
        .taxa
        .splitCsv(header: true)
        .map { meta, data -> 
            [meta, data["predicted.name"].replace(" ", "_")] 
        }
        .join(
            ch_manifest
                .filter { !it.meta.taxa }
                .map { [it.meta, it] }
        )
        .map { meta, taxa, data -> 
            data.meta.taxa = taxa
            return data 
        }
        .concat(ch_manifest.filter { it.meta.taxa })

    /* 
    =============================================================================================================================
        CLUSTER ASSIGNMENT (FLOC)
    =============================================================================================================================
    */

    // Group samples by taxa for clustering
    // FIXED: Check actual database state instead of relying on cached metadata
    ch_taxa_groups = ch_manifest_with_taxa
        .map { data ->
            // Check if this sample's signature already exists in the database
            def sig_file = file(params.db) / data.meta.taxa / 'sig' / "${data.meta.id}.sig.gz"
            def already_processed = sig_file.exists()
            
            [data, already_processed]
        }
        .filter { data, already_processed -> 
            !already_processed  // Only keep samples that haven't been processed
        }
        .map { data, already_processed -> data }  // Remove the flag, keep the data
        .map { [
            [taxa: it.meta.taxa], 
            it.assembly,
            it.meta.timestamp
        ] }
        .groupTuple(by: 0, sort: true)
        .map { meta, assemblies, timestamps -> 
            meta.timestamp = timestamps.min()
            def sig_db = file(params.db) / meta.taxa / 'sig'
            [meta, assemblies, sig_db.exists() ? sig_db : []]
        }

    // Run FLOC clustering
    FLOC(ch_taxa_groups)

    // Parse cluster assignments
    ch_clusters = FLOC
        .out
        .clusters
        .splitCsv(header: true)
        .transpose()
        .map { meta, csv -> 
            [csv.sid.toString(), csv.cluster] 
        }

    // Enrich metadata with cluster assignments
    ch_manifest_final = ch_manifest_with_taxa
        .filter { !it.meta.cluster }
        .map { [it.meta.id.toString(), it] }
        .join(ch_clusters, by: 0)
        .map { id, data, cluster ->
            data.meta.cluster = cluster
            return data
        }
        .concat(ch_manifest_with_taxa.filter { it.meta.cluster })

    // Prepare cluster-level distance matrices
    ch_db_dist = FLOC
        .out
        .dist
        .transpose()
        .map { meta, result -> 
            def cluster_name = result.name.replace('.tsv.gz', '')
            [
                [taxa: meta.taxa, cluster: cluster_name, timestamp: meta.timestamp], 
                result
            ]
        }

    emit:
    db       = FLOC.out.sigs      // channel: [ val(meta), path(sigs) ]
    db_dist  = ch_db_dist          // channel: [ val(meta), path(dist) ]
    manifest = ch_manifest_final   // channel: [ val(meta), path(assembly), ... ]
    versions = ch_versions         // channel: [ path(versions.yml) ]
}
