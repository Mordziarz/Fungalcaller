process SUMMARIZE_BUSCO {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path(busco_dirs)

    output:
    path("busco_summary_table.csv"), emit: csv

    script:
    """
    python3 ${workflow.projectDir}/bin/summarize_busco.py .
    """
}
