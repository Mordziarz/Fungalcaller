#!/usr/bin/env python3
import os
import sys
import subprocess
import csv
import glob


def load_known_samples():
    """
    Loads sample IDs from samples_table.csv if present in common candidate locations.
    Returns a list sorted by descending length so longer sample names match first.
    """
    candidates = [
        "samples_table.csv",
        os.path.join(os.getcwd(), "samples_table.csv"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "samples_table.csv"),
    ]
    candidates.extend(glob.glob(os.path.join(os.getcwd(), "..", "..", "*", "samples_table.csv")))

    samples = set()
    for path in candidates:
        abs_path = os.path.abspath(path)
        if os.path.isfile(abs_path):
            try:
                with open(abs_path, "r", newline="", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    if reader.fieldnames and "sample" in reader.fieldnames:
                        for row in reader:
                            val = (row.get("sample") or "").strip()
                            if val:
                                samples.add(val)
                if samples:
                    break
            except OSError:
                continue

    return sorted(samples, key=len, reverse=True)


def extract_sample_id(header, known_samples):
    """
    Extracts a clean taxon/sample identifier from a FASTA header without breaking
    underscore-delimited IDs such as GCA_000240135.3_ASM24013v3_genomic.
    """
    token = header.strip().split()[0] if header.strip() else ""
    for sample in known_samples:
        if token == sample or token.startswith(sample + "_") or token.startswith(sample + "|") or token.startswith(sample + ":"):
            return sample
    for sample in known_samples:
        if sample in token:
            return sample
    return token


def parse_fasta_seqs(input_fa, known_samples):
    """
    Parses a FASTA file into a dict mapping sample_id -> aligned sequence string.
    """
    seqs = {}
    curr_sample = None
    curr_seq = []

    with open(input_fa, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_sample is not None and curr_seq:
                    seqs[curr_sample] = "".join(curr_seq)
                curr_sample = extract_sample_id(line[1:], known_samples)
                curr_seq = []
            else:
                curr_seq.append(line)

        if curr_sample is not None and curr_seq:
            seqs[curr_sample] = "".join(curr_seq)

    return seqs


def analyze_gene_alignment(gene_name, seqs, standard_samples):
    """
    Analyzes a single gene alignment for length, mutations (variability), and gaps (indels).
    """
    max_len = max((len(s) for s in seqs.values()), default=0)
    ordered_seqs = [seqs.get(s, "").ljust(max_len, "-") for s in standard_samples]

    num_taxa = len(standard_samples)
    aln_length = max_len

    variable_sites = 0
    parsimony_informative = 0
    gap_columns = 0
    total_gaps = 0

    for col_idx in range(aln_length):
        col_chars = [seq[col_idx] for seq in ordered_seqs]

        if "-" in col_chars:
            gap_columns += 1
            total_gaps += col_chars.count("-")

        non_gap_chars = [c for c in col_chars if c != "-"]

        if len(set(non_gap_chars)) > 1:
            variable_sites += 1
            char_counts = {c: non_gap_chars.count(c) for c in set(non_gap_chars)}
            multiton_count = sum(1 for count in char_counts.values() if count >= 2)
            if multiton_count >= 2:
                parsimony_informative += 1

    return {
        "gene": gene_name,
        "taxa_count": num_taxa,
        "alignment_length": aln_length,
        "variable_sites": variable_sites,
        "parsimony_informative_sites": parsimony_informative,
        "gap_columns": gap_columns,
        "total_gaps": total_gaps,
        "fraction_gaps": round(total_gaps / (num_taxa * aln_length), 4) if (num_taxa * aln_length) > 0 else 0,
    }


def write_normalized_fasta(seqs, output_fa, standard_samples):
    """
    Writes a normalized FASTA containing all standard_samples padded to max_len with '-'.
    """
    max_len = max((len(s) for s in seqs.values()), default=0)
    with open(output_fa, "w", encoding="utf-8") as out:
        for sample in standard_samples:
            seq = seqs.get(sample, "")
            if max_len > 0:
                padded = seq.ljust(max_len, "-")[:max_len]
                out.write(f">{sample}\n{padded}\n")


def main():
    input_files = sys.argv[1:]
    if not input_files:
        print("[ERROR] No input files provided for PhyKIT!", file=sys.stderr)
        sys.exit(1)

    known_samples = load_known_samples()
    if known_samples:
        print(f"[INFO] Loaded {len(known_samples)} sample names from samples_table.csv")

    parsed_alignments = []
    all_samples = set(known_samples)

    for fa in input_files:
        seqs = parse_fasta_seqs(fa, known_samples)
        parsed_alignments.append((fa, seqs))
        all_samples.update(seqs.keys())

    standard_samples = sorted(all_samples)
    print(f"[INFO] Established fixed sample sorting order ({len(standard_samples)} taxa): {standard_samples}")

    cleaned_dir = "cleaned_for_phykit"
    os.makedirs(cleaned_dir, exist_ok=True)

    gene_stats_list = []
    cleaned_list_file = "alignment_list_cleaned.txt"

    with open(cleaned_list_file, "w", encoding="utf-8") as list_out:
        for fa, seqs in parsed_alignments:
            base_name = os.path.basename(fa)
            gene_name = base_name.replace("_trimmed.fa", "").replace(".fa", "")
            out_fa = os.path.join(cleaned_dir, base_name)

            write_normalized_fasta(seqs, out_fa, standard_samples)
            stats = analyze_gene_alignment(gene_name, seqs, standard_samples)
            gene_stats_list.append(stats)

            list_out.write(os.path.abspath(out_fa) + "\n")

    report_file = "gene_variation_report.csv"
    with open(report_file, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "gene",
            "taxa_count",
            "alignment_length",
            "variable_sites",
            "parsimony_informative_sites",
            "gap_columns",
            "total_gaps",
            "fraction_gaps",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in gene_stats_list:
            writer.writerow(row)
    print(f"[SUCCESS] Generated gene variation report: {report_file}")

    cmd = ["phykit", "create_concat", "-a", cleaned_list_file, "-p", "concatenated_output"]
    print(f"[INFO] Running PhyKIT: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] PhyKIT encountered an error: {e}", file=sys.stderr)
        sys.exit(1)

    if os.path.exists("concatenated_output.fa"):
        os.replace("concatenated_output.fa", "concatenated_alignment.fasta")
    if os.path.exists("concatenated_output.partition"):
        os.replace("concatenated_output.partition", "partitions.txt")


if __name__ == "__main__":
    main()
