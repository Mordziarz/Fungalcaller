process JELLYFISH {
    tag "$sample_id"
    publishDir "${params.outdir}/jellyfish", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_histo.txt"), emit: histo

    script:
    """
    if [[ "${fasta}" == *.gz ]]; then
        zcat ${fasta} > input_seqs.fa
    else
        cp ${fasta} input_seqs.fa
    fi

    jellyfish count -m 21 -s 100M -t ${task.cpus} -C input_seqs.fa -o ${sample_id}_mer_counts.jf

    jellyfish histo ${sample_id}_mer_counts.jf > ${sample_id}_histo.txt

    rm input_seqs.fa ${sample_id}_mer_counts.jf
    """
}
