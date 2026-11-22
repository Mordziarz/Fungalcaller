process ALIGNMENT_CONCATENATION {
    tag "concatenation"
    label 'process_medium'
    publishDir "${params.outdir}/concatenated_alignment", mode: 'copy'

    input:
    path alignments
    path samples_csv

    output:
    path "concatenated_alignment.fasta", emit: alignment
    path "partitions.txt", emit: partitions
    path "gene_variation_report.csv", emit: report

    script:
    """
    set -euo pipefail

    echo "Concatenating normalized orthologs and auditing variants..."

    python3 ${workflow.projectDir}/bin/run_phykit_concat.py --samples-table "${samples_csv}" ${alignments}
    """
}
