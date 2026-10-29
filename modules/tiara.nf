process TIARA {
    tag "$sample_id"
    publishDir "${params.outdir}/tiara", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_clean.fasta"), emit: cleaned_fasta
    path("${sample_id}_tiara.txt"), emit: report

    script:
    """
    tr -d '\\r' < ${fasta} > input_clean_${sample_id}.fasta

    tiara -i input_clean_${sample_id}.fasta \
          -o ${sample_id}_tiara.txt \
          -t ${task.cpus} \
          --to_fasta euk mit unk \
          -m 1000 \
          -p 0.5

    final_fasta="${sample_id}_clean.fasta"
    > "\$final_fasta"

    for cat_file in eukarya_* mitochondrion_* unknown_*; do
        if [ -f "\$cat_file" ]; then
            echo "Adding fraction: \$cat_file"
            cat "\$cat_file" >> "\$final_fasta"
        fi
    done

    rm -f bacteria_* archaea_* prokarya_*

    if [ -f "\$final_fasta" ]; then
        ORIG_COUNT=\$(grep -c "^>" ${fasta} || true)
        CLEAN_COUNT=\$(grep -c "^>" "\$final_fasta" || true)
        echo "Original sequences: \$ORIG_COUNT"
        echo "After Tiara cleaning (Euk+Mit+Unk >= 1000bp): \$CLEAN_COUNT"
    else
        echo "Error: Output file \$final_fasta was not created!"
        exit 1
    fi
    """
}
