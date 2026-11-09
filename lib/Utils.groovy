class Utils {

    static Map cluster_meta(rec) {
        def CM_KEYS = ['taxa', 'cluster', 'timestamp']

        if (!rec) {
            error "cluster_meta: record is empty or null"
        }

        if (!(rec instanceof Map)) {
            error "cluster_meta: record must be map"
        }

        // Fail if required CM keys are missing
        def missing = CM_KEYS.findAll { !rec.containsKey(it) }
        if (missing) {
            error "cluster_meta: missing required CM keys: ${missing}"
        }

        def out = new LinkedHashMap()

        // Insert CM keys in defined order
        CM_KEYS.each { k ->
            out[k] = rec[k]
        }

        return out
    }

    // Make sure an existing taxa database was built with the same cluster distance as this run.
    // Fails if config.json is missing or clust_dist differs. A different clust_assign_hashfreq only warns.
    static void check_db_config(config_file, String taxa, Map params, log) {
        def db_dir = config_file.parent

        if (!config_file.exists()) {
            nextflow.Nextflow.error(
                "The existing ${taxa} database (${db_dir}) has no config.json, so the cluster distance it was built with is unknown. " +
                "Rebuilding the database with the current defaults is strongly recommended (see the v2.2.0 CHANGELOG). " +
                "To keep using it as it is, create ${config_file} containing {\"clust_dist\": <value>} with the distance it was built with " +
                "(0.03 for BigBacter v2.0.x-v2.1.x defaults) and run again."
            )
        }

        def config
        try {
            config = new groovy.json.JsonSlurper().parseText(config_file.text)
        } catch (Exception e) {
            nextflow.Nextflow.error("Could not read ${config_file}: ${e.message}")
        }

        if (config?.clust_dist == null) {
            nextflow.Nextflow.error("${config_file} does not record clust_dist.")
        }

        if (!same_number(config.clust_dist, params.clust_dist)) {
            nextflow.Nextflow.error(
                "--clust_dist ${params.clust_dist} does not match the distance the ${taxa} database was built with " +
                "(${config.clust_dist}, from ${config_file}). Run again with --clust_dist ${config.clust_dist}, or use a new database."
            )
        }

        if (config.clust_assign_hashfreq != null && !same_number(config.clust_assign_hashfreq, params.clust_assign_hashfreq)) {
            log.warn "--clust_assign_hashfreq ${params.clust_assign_hashfreq} differs from the value the ${taxa} database was built with " +
                     "(${config.clust_assign_hashfreq}). Cluster assignments may be less consistent with earlier runs."
        }
    }

    // Compare two numbers given as numbers or strings, so 0.1, 0.10 and 1e-1 are equal
    static boolean same_number(a, b) {
        try {
            return new BigDecimal(a.toString()).compareTo(new BigDecimal(b.toString())) == 0
        } catch (NumberFormatException e) {
            nextflow.Nextflow.error("Expected a number but got '${a}' and '${b}'.")
        }
    }

}

