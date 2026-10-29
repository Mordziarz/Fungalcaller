process ANTISMASH {
    tag "$sample_id"
    publishDir "${params.outdir}/antismash", mode: 'copy'

    input:
    tuple val(sample_id), path(genbank)

    output:
    tuple val(sample_id), path("${sample_id}_antismash_out"), emit: results

    script:
    """
    set -euo pipefail

    export HOME="\$PWD/tmp_home"
    export MPLCONFIGDIR="\$PWD/matplotlib_config"

    mkdir -p "\$HOME" "\$MPLCONFIGDIR"

    echo "=========================================="
    echo "antiSMASH: sample ${sample_id}"
    echo "Input file: ${genbank}"
    echo "=========================================="

    ANTISMASH_DB_DIR="/opt/conda/envs/fungi/lib/python3.9/site-packages/antismash/databases"

    if [ -d "/databases/antismash" ]; then
        for db in /databases/antismash/*; do
            db_name=\$(basename "\$db")
            if [ ! -e "\${ANTISMASH_DB_DIR}/\${db_name}" ]; then
                ln -s "\$db" "\${ANTISMASH_DB_DIR}/\${db_name}" || true
            fi
        done
    fi

    python3 - "${genbank}" annotations_for_antismash.gbk <<'PY'
    import sys
    from Bio import SeqIO

    input_gbk = sys.argv[1]
    output_gbk = sys.argv[2]

    records = list(SeqIO.parse(input_gbk, "genbank"))

    if not records:
        sys.exit(f"ERROR: File {input_gbk} does not contain GenBank records.")

    valid_records = []
    gene_counter = 1

    for record in records:
        has_cds = False
        for feature in record.features:
            if feature.type in ["CDS", "gene", "mRNA"]:
                has_cds = True
                unique_name = f"g_{gene_counter:06d}"
                
                feature.qualifiers["locus_tag"] = [unique_name]
                if "protein_id" in feature.qualifiers:
                    feature.qualifiers["protein_id"] = [unique_name]
                if "gene" in feature.qualifiers:
                    feature.qualifiers["gene"] = [unique_name]
                
            if feature.type == "CDS":
                gene_counter += 1
                
        if has_cds:
            valid_records.append(record)

    if not valid_records:
        sys.exit("ERROR: No GenBank records contain CDS annotations.")

    SeqIO.write(valid_records, output_gbk, "genbank")
    print(f"Saved {len(valid_records)} records with unique gene IDs to {output_gbk}")
    PY

    antismash \
        --taxon fungi \
        --cb-general \
        --cb-knownclusters \
        --asf \
        --output-dir "${sample_id}_antismash_out" \
        annotations_for_antismash.gbk

    if [ ! -d "${sample_id}_antismash_out" ]; then
        echo "ERROR: antiSMASH did not create the output directory."
        exit 1
    fi

    echo "antiSMASH successfully completed analysis for sample ${sample_id}."
    """
}
