process BUSCO {
    tag "$sample_id"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path("${sample_id}_busco_out"), emit: report

    script:
    """
    busco -i ${fasta} \
          -m genome \
          -l ${params.busco_lineage} \
          -o ${sample_id}_busco_out \
          -c ${task.cpus} \
          --quiet -f
    """
}
