process TRIMAL {
    tag "${aligned_fa.baseName}"
    label 'process_low'
    container 'mordziarz/fungi_pipeline:latest'
    
    publishDir "${params.outdir}/mafft_single_copy", mode: 'copy'

    input:
    path(aligned_fa)

    output:
    path("*_trimmed.fa"), emit: trimmed_fa

    script:
    """
    set -euo pipefail
    trimal -in ${aligned_fa} -out ${aligned_fa.baseName}_trimmed.fa -automated1
    """
}
