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

}

