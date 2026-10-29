process INTERPROSCAN {
    tag "$id"
    container 'mordziarz/fungi_pipeline:latest'

    publishDir "${params.outdir}/interproscan", mode: 'copy'

    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("${id}.tsv"), emit: tsv

    script:
    """
    set -euo pipefail

    python3 - <<'PY'
    from Bio import SeqIO
    from Bio.Seq import Seq

    input_file = "${fasta}"
    cleaned_file = "clean_${id}.faa"

    valid_aas = set("ACDEFGHIKLMNPQRSTVWYXUOBZJ")

    records = []
    for rec in SeqIO.parse(input_file, "fasta"):
        rec.id = rec.id.replace(".", "_").replace(" ", "_")
        rec.description = ""

        seq_str = str(rec.seq).upper()
        cleaned_seq = "".join([char if char in valid_aas else 'X' for char in seq_str])

        if len(cleaned_seq) >= 5:
            rec.seq = Seq(cleaned_seq)
            records.append(rec)

    SeqIO.write(records, cleaned_file, "fasta")
    print(f"Cleaned FASTA file for InterProScan: {len(records)} sequences.")
    PY

    interproscan.sh -i clean_${id}.faa \
                    -f tsv \
                    -dp \
                    -appl Pfam,SuperFamily,Gene3D,SMART,ProSiteProfiles,ProSitePatterns \
                    -cpu ${task.cpus} \
                    -b ${id}
    """
}
