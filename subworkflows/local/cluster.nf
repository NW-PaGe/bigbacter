//
// Assign PopPUNK clusters for each isolate
//

// Modules
include { BB_CLUSTER     } from '../../modules/local/bb-cluster'
include { POPPUNK_ASSIGN } from '../../modules/local/poppunk-assign'
include { POPPUNK_VISUAL } from '../../modules/local/poppunk-visualize'
include { RESOLVE_MERGED } from '../../modules/local/resolve-merged-clusters'

/*
=============================================================================================================================
    SUBWORKFLOW
=============================================================================================================================
*/
workflow CLUSTER {
    take:
    ch_man        // channel: [ val(sample), val(taxa), file(assembly) ]
    ch_timestamp  // channel: val(timestamp)

    main:
    ch_versions = Channel.empty()
    /*
    =============================================================================================================================
        ASSIGN CLUSTERS
        - Determine if a database exists for each species (PopPUNK or bb-cluster; default: bb-cluster)
        - Return most recent database version (if a database exists)
        - Assign clusters with selected database & method
    =============================================================================================================================
    */
    // Fetch current database info
    ch_man
        .groupTuple( by: 1 )
        .combine( ch_timestamp )
        .map{ sample, taxa, assembly, timestamp -> [ taxa: taxa, sample: sample, assembly: assembly ] + getLastDb( taxa, timestamp ) }
        .branch{ it ->
            pp: it.db_type == 'pp_db'
            bb: it.db_type == 'bb_db'
        }
        .set { ch_db }

    /*
    ====================
        PopPUNK
    ====================
    */
    // MODULE: Assign clusters with PopPUNK
    POPPUNK_ASSIGN (
        ch_db.pp.map{ [ it.taxa, it.sample, it.assembly, it.db_file ] },
        ch_timestamp
    )
    ch_versions = ch_versions.mix(POPPUNK_ASSIGN.out.versions)

    // MODULE: Create visuals for new PopPUNK database
    POPPUNK_VISUAL(
        POPPUNK_ASSIGN.out.db,
        ch_timestamp
    )
    ch_versions = ch_versions.mix(POPPUNK_VISUAL.out.versions)

    /*
    ====================
        bb-cluster
    ====================
    */
    // MODULE: Assign clusters with bb-cluster
    BB_CLUSTER (
        ch_db.bb.map{ [ it.taxa, it.assembly, it.db_file ] },
        ch_timestamp
    )
    ch_versions = ch_versions.mix(BB_CLUSTER.out.versions)

    // Load cluster results
    POPPUNK_ASSIGN
        .out
        .cluster_results
        .concat( BB_CLUSTER.out.clusters )
        .splitCsv( header: true, elem: 1 )
        .transpose()
        .map{ taxa, data -> [ taxa, 
                              data.containsKey('name') ? data['name'].tokenize('.')[0] : data['Taxon'], 
                              data.containsKey('cluster') ? data['cluster'] : data['Cluster'] ] }
        .join( ch_man.map { sample, taxa, assembly  -> [ taxa, sample ] }, by: [0,1] )
        .map{ taxa, sample, cluster -> [ taxa: taxa, sample: sample, cluster: cluster ] }
        .set{ ch_clusters }

    // Create channnel for new database files
    BB_CLUSTER
        .out
        .db
        .map{ taxa, db -> [ taxa, db, 'bb_db' ] }
        .concat(POPPUNK_ASSIGN.out.db.map{ taxa, db -> [ taxa, db, 'pp_db' ] })
        .set{ ch_new_db }

    // Create channel for pairwise distances used for clustering
    BB_CLUSTER
        .out
        .dist
        .concat( POPPUNK_ASSIGN.out.core_acc_dist )
        .set{ ch_dist }
    
    /*
    =============================================================================================================================
        RESOLVE MERGED CLUSTERS
        - Only applies to PopPUNK (all samples still passed through)
        - Merged clusters are denoted by an underscore ('_') in the cluster name
    =============================================================================================================================
    */
    // MODULE: Resolve merged clusters
    RESOLVE_MERGED (
        ch_clusters
            .filter{ it.cluster.contains('_') }
            .map{ [ it.taxa, it.sample ] }
            .combine( POPPUNK_ASSIGN.out.jaccard_dist, by: 0 )
    )

    RESOLVE_MERGED
        .out
        .neighbor
        .map{ taxa, sample, neighbor -> [ taxa, neighbor ] }
        .distinct()
        .groupTuple(by: 0)
        .map{ resolveMerged( it ) }
        .combine( RESOLVE_MERGED.out.neighbor.map{ taxa, sample, neighbor -> [ taxa, neighbor, sample ] }, by: [0,1] )
        .map{ taxa, neighbor, cluster, sample -> [ taxa, sample, cluster ] }
        .set{ ch_resolved }
    ch_clusters
        .map{ [ it.taxa, it.sample, it ] }
        .join( ch_resolved, by: [0,1], remainder: true)
        .map{ taxa, sample, data, resolved -> data.cluster = resolved? resolved : data.cluster
                                              data }
        .set{ ch_clusters }

    /*
    =============================================================================================================================
        FORMAT CLUSTERS
    =============================================================================================================================
    */
    // Add padding to cluster numbers
    ch_clusters
        .map{ it.cluster = it.cluster.padLeft(5, "0")
              it }
        .set{ ch_clusters }

    emit:
    clusters = ch_clusters // channel: [ val(sample), val(cluster) ]
    new_db   = ch_new_db   // channel: [ val(taxa), path(new_pp_db)]
    dist     = ch_dist     // channel: [ val(taxa), path(dist) ]
    versions = ch_versions // channel: [ versions.yml ]
}

/*
=============================================================================================================================
    FUNCTIONS
=============================================================================================================================
*/
// Function for determining the most recent database
def getLastDb ( taxa, timestamp ) {
    def taxon_path = file(params.db, checkIfExists: false).resolve(taxa)
    def db_opts = ['pp_db', 'bb_db']
    def db_new  = true
    def db_type = 'bb_db'

    if( taxon_path.exists() ){
        db_type = taxon_path.list().findAll{ it in db_opts }
        if( db_type.size() > 1 ){
            msg = """
            Databases can only use one clustering method. Please remove one of the following directories to proceed:

            ${ db_type.collect{ taxon_path.resolve( it ) }.join('\n') }

            """
            exit 1, msg
        }

        db_path   = taxon_path.resolve( db_type[0] )
        db_files  = db_path.list()
        db_new    = db_files.size() == 0 ? true : false

    }
    db_type   = db_new ? 'bb_db' : db_type[0]
    db_path   = taxon_path.resolve( db_type )
    last_file = db_new ? [] : db_path.resolve( db_files.sort().last() )
        
    return [ db_type: db_type, db_path: db_path, db_file: last_file ]
}

// Function for resolving merged PopPUNK clusters
def resolveMerged( data ){
    def taxa      = data[0]
    def neighbors = data[1]
    def result = file(params.db)
        .resolve(taxa)
        .resolve('clusters')
        .list()
        .collectMany{ c -> 
            file( params.db )
                .resolve( taxa )
                .resolve( 'clusters' )
                .resolve( c )
                .resolve('snippy')
                .list()
                .findAll{ f -> 
                    f.replaceAll('.tar.gz', '').contains( neighbors ) 
                }.collect{ f -> [ taxa, f.replaceAll('.tar.gz', ''), c ] }
        }
    return result
}