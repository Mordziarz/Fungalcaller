process GFF_TO_GBK {
    tag "$sample_id"
    publishDir "${params.outdir}/genbank", mode: 'copy'

    input:
    tuple val(sample_id), path(gff), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.gbk"), emit: gbk

    script:
    """
    python3 - "${fasta}" "${gff}" "${sample_id}.gbk" <<'PY'
    import sys
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    fasta_file = sys.argv[1]
    gff_file = sys.argv[2]
    output_gbk = sys.argv[3]

    records = {rec.id: rec for rec in SeqIO.parse(fasta_file, "fasta")}
    for rec in records.values():
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "linear"

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attrs = parts
            if feature_type != "CDS":
                continue
            if seqid not in records:
                continue

            start_idx = int(start) - 1
            end_idx = int(end)
            strand_val = 1 if strand == "+" else -1

            gene_id = "cds"
            for attr in attrs.split(";"):
                if attr.startswith("gene_id=") or attr.startswith("ID=") or attr.startswith("Parent="):
                    gene_id = attr.split("=")[1]
                    break

            cds_feat = SeqFeature(
                FeatureLocation(start_idx, end_idx, strand=strand_val),
                type="CDS",
                qualifiers={"locus_tag": [gene_id], "transl_table": [1]}
            )
            records[seqid].features.append(cds_feat)

    SeqIO.write(list(records.values()), output_gbk, "genbank")
    print(f"Generated GenBank file: {output_gbk}")
    PY
    """
}
