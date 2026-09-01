process GAMBIT_STAGE {

    input:
    val ids
    path gambit_gdb, stageAs: "db/*"
    path gambit_gs,  stageAs: "db/*"

    output:
    path "db/", emit: db

    script:
    """
    sleep 1
    """
}