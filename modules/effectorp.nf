process EFFECTORP {
    tag "$id"
    label 'process_medium'

    container 'mordziarz/fungi_pipeline:latest'

    publishDir "${params.outdir}/effectorp", mode: 'copy'

    input:
    tuple val(id), path(protein_fasta)

    output:
    tuple val(id), path("${id}_effectorp_report.txt"), emit: report

    script:
    """
    echo "Running EffectorP 3.0 for sample: $id on the whole proteome..."

    python3 /opt/EffectorP-3.0/EffectorP.py -i ${protein_fasta} > ${id}_effectorp_report.txt

    if [ ! -s "${id}_effectorp_report.txt" ]; then
        echo "ERROR: EffectorP report is empty."
        exit 1
    fi

    echo "EffectorP 3.0 finished successfully for sample: $id"
    """
}
