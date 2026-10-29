process SPADES {
    tag "$sample_id"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta"), emit: contigs

    script:
    """
    spades.py -1 ${r1} -2 ${r2} -o spades_out --careful
    cp spades_out/contigs.fasta ${sample_id}_contigs.fasta
    """
}
