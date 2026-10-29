process RAXML {
    tag "raxml_phylogeny"
    label 'process_high'
    container 'mordziarz/fungi_pipeline:latest'
    publishDir "${params.outdir}/raxml", mode: 'copy'

    input:
    path(alignment)

    output:
    path("raxml_out/*"), emit: results

    script:
    """
    set -euo pipefail
    mkdir -p raxml_out

    echo "Running RAxML-NG using ${task.cpus} cores on file: ${alignment}..."

    raxml-ng --all \
        --msa ${alignment} \
        --model LG+G4 \
        --prefix raxml_out/species_tree \
        --threads ${task.cpus} \
        --bs-trees 1000 \
        --outgroup ${params.outgroup}
    """
}
