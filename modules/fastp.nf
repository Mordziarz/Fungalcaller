process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_val_1.fastq.gz"), path("${sample_id}_val_2.fastq.gz"), emit: reads
    path("${sample_id}_fastp.json")
    path("${sample_id}_fastp.html")

    script:
    """
    fastp -i ${r1} -I ${r2} \
          -o ${sample_id}_val_1.fastq.gz -O ${sample_id}_val_2.fastq.gz \
          --json ${sample_id}_fastp.json --html ${sample_id}_fastp.html
    """
}
