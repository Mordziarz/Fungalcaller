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

    orthofinder -f ortho_input -t ${task.cpus}

    mv ortho_input/OrthoFinder/Results_* orthofinder
    """
}
