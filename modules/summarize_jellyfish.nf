process SUMMARIZE_JELLYFISH {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path(histo_files) // Zbiera wszystkie pliki .histo ze wszystkich próbek

    output:
    path("jellyfish_summary_table.csv"), emit: csv

    script:
    """
    python3 ${workflow.projectDir}/bin/summarize_jellyfish.py .
    """
}
