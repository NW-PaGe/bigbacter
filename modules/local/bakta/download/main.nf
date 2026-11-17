process BAKTA_DOWNLOAD {
    label "process_low"

    input:
    val type // 'light' or 'full'

    output:
    path "db-${type}",   emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # downloads the newest database the installed Bakta version can use, and sets up AMRFinderPlus
    bakta_db download --output download --type ${type} ${args}

    # bakta_db names the directory 'db-light' or 'db'
    mv download/db* db-${type}

    # replace symlinks (AMRFinderPlus 'latest') with what they point to, as object stores cannot hold them
    find db-${type} -type l | while read -r link; do
        target=\$(readlink -f "\$link")
        rm "\$link"
        cp -r "\$target" "\$link"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: "\$(bakta --version 2>/dev/null | sed -n 's/^bakta //p' | tail -n 1)"
    END_VERSIONS
    """

    stub:
    """
    mkdir db-${type}
    echo '{"major": 0, "minor": 0, "type": "${type}"}' > db-${type}/version.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: stub
    END_VERSIONS
    """
}
