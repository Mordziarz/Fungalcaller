process REPEATMASKER {
    tag "$sample_id"
    publishDir "${params.outdir}/repeatmasker", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${fasta}.masked"), emit: masked_fasta
    path("*.gff"), emit: gff
    path("*.tbl"), emit: tbl
    path("*.out"), emit: out

    script:
    """
    RepeatMasker -species fungi \
                 -xsmall \
                 -gff \
                 -pa ${task.cpus} \
                 ${fasta}
    """
}
