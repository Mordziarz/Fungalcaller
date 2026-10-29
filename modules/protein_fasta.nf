process PROTEIN_FASTA {
    tag "$sample_id"
    publishDir "${params.outdir}/sequences", mode: 'copy'

    input:
    tuple val(sample_id), path(gff), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.faa"), emit: proteins
    tuple val(sample_id), path("${sample_id}_cds.fasta"), emit: cds

    script:
    """
    gffread ${gff} -g ${fasta} -y ${sample_id}_raw.faa

    gffread ${gff} -g ${fasta} -x ${sample_id}_cds_raw.fasta

    awk -v p="${sample_id}" '
        /^>/ { 
            \$0 = ">" p "_" substr(\$0, 2) 
        } 
        !/^>/ { 
            gsub(/[.]/, "", \$0) 
        } 
        {print}
    ' ${sample_id}_raw.faa > ${sample_id}.faa

    awk -v p="${sample_id}" '/^>/ {\$0 = ">" p "_" substr(\$0, 2)} {print}' ${sample_id}_cds_raw.fasta > ${sample_id}_cds.fasta
    """
}
