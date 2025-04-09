process POPPUNK_ASSIGN {
    tag "${taxa}"
    label 'process_high'

    input:
    tuple val(taxa), val(sample), path(assembly), path(db)
    val timestamp

    output:
    tuple val(taxa), path('clusters.csv'),                  emit: cluster_results
    tuple val(taxa), path('pp-core-acc-dist.txt.gz'),       emit: core_acc_dist
    tuple val(taxa), path('pp-jaccard-dist.txt.gz'),        emit: jaccard_dist
    tuple val(taxa), path('*.tar.gz', includeInputs: true), emit: db
    path 'versions.yml',                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args           = task.ext.args ?: ''
    prefix         = "${timestamp}-${taxa}"
    db_name        = db.name
    '''
    # decompress database
    db_comp="!{db_name}"
    db=${db_comp%.tar.gz}
    tar -xzvf !{db} -C ./

    #### CREATE QFILE ####
    echo !{sample.join(',')} | tr ',' '\n' > s_col
    echo !{assembly.join(',')} | tr ',' '\n' > a_col
    paste s_col a_col > qfile.txt

    #### ASSIGN CLUSTERS ####
    # run PopPUNK
    echo -e "\nAssigning clusters with PopPUNK.\n"
    poppunk_assign \
        !{args} \
        --db ${db}  \
        --query qfile.txt \
        --output !{timestamp} \
        --threads !{task.cpus} \
        --max-zero-dist !{params.pp_maxzerodist} \
        --max-merge !{params.pp_maxmerge} \
        --max-a-dist !{params.pp_maxadist} \
        --max-pi-dist !{params.pp_maxpidist} \
        !{params.pp_runqc ? '--run-qc' : ''}


    # check that all samples made it into the new database
    for SAMPLE in $(cat qfile.txt | cut -f 1)
    do
        STATUS=$(cat !{timestamp}/!{timestamp}_clusters.csv | awk -v FS=',' -v sample="${SAMPLE}" '$1==sample {print "true"}')
        if [[ $STATUS != "true" ]]
        then
            echo "ERROR: $SAMPLE was not found in !{timestamp}/!{timestamp}_clusters.csv. This is likely due to the sample failing the PopPUNK QC. Try running BigBacter with '--pp_runqc false' to bypass this." && exit 1
        fi
    done

    # compress new database (output) & remove old database
    tar -czvf !{timestamp}.tar.gz !{timestamp} && rm -r ${db}*

    #### CALCULATE DISTANCES ####
    # core & accessory distances
    sketchlib query dist --cpus !{task.cpus} */*[^refs].h5 | gzip > pp-core-acc-dist.txt.gz
    sketchlib query jaccard --cpus !{task.cpus} */*[^refs].h5 | gzip > pp-jaccard-dist.txt.gz

    #### RENAME CLUSTER RESULTS ####
    cp */*[^_unword]_clusters.csv clusters.csv
    
    # version info
    echo "!{task.process}:\n    poppunk: $(poppunk_assign --version | cut -f 2 -d ' ' | tr -d '\n\r\t ')" > versions.yml
    '''
}
