process ORTHOFINDER {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(proteins)

    output:
    path("orthofinder"), emit: results

    script:
    """
    mkdir -p ortho_input
    cp ${proteins} ortho_input/

    A_THREADS=${task.cpus}
    [ "\$A_THREADS" -gt 8 ] && A_THREADS=8
    [ "\$A_THREADS" -lt 1 ] && A_THREADS=1

    T_THREADS=\$(expr ${task.cpus} / \$A_THREADS)
    [ "\$T_THREADS" -lt 1 ] && T_THREADS=1

    orthofinder -f ortho_input -t "\$T_THREADS" -a "\$A_THREADS"

    mv ortho_input/OrthoFinder/Results_* orthofinder
    """
}
