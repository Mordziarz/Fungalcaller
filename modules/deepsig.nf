process DEEPSIG {
    tag "$sample_id"
    publishDir "${params.outdir}/deepsig", mode: 'copy'

    input:
    tuple val(sample_id), path(protein_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_deepsig.txt"), emit: results

    script:
    """
    set -euo pipefail

    echo "Running local DeepSig for sample: ${sample_id}"

    deepsig -f "${protein_fasta}" -o "${sample_id}_deepsig.txt" -k euk

    if [ ! -s "${sample_id}_deepsig.txt" ]; then
        echo "ERROR: Local DeepSig did not generate a valid output file."
        exit 1
    fi

    echo "DeepSig successfully finished running locally for sample: ${sample_id}"
    """
}
