process AUGUSTUS {
    tag "$sample_id"
    publishDir "${params.outdir}/augustus", mode: 'copy'

    input:
    tuple val(sample_id), path(masked_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_augustus.gff3"), path("clean_masked.fasta"), emit: gff_and_fasta

    script:
    """
    # 1. Czyszczenie znaków Windows
    tr -d '\\r' < ${masked_fasta} > clean_masked.fasta

    # 2. Uruchomienie Augustusa
    augustus --species=fusarium \
             --gff3=on \
             --softmasking=on \
             clean_masked.fasta > ${sample_id}_augustus.gff3
    """
}
