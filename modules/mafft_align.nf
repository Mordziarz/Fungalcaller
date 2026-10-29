process MAFFT_ALIGN {
    tag "${fasta.baseName}"
    container 'mordziarz/fungi_pipeline:latest'
    input:
    path(fasta)

    output:
    path("${fasta.baseName}_aligned.fa"), emit: aligned_fa

    script:
    """
    set -euo pipefail
    mafft --auto --thread ${task.cpus} ${fasta} > ${fasta.baseName}_aligned.fa
    """
}
