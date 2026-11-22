#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import sys

NON_INFORMATIVE = {"-", "N", "?", "X"}


def load_known_samples(samples_csv_path):
    """
    Loads exact sample IDs from the provided CSV file (column 'sample').
    Returns a list sorted by descending length so longer sample names match first.
    """
    if not samples_csv_path or not os.path.isfile(samples_csv_path):
        print(f"[WARN] Samples table not found at '{samples_csv_path}'. Falling back to FASTA headers.", file=sys.stderr)
        return []

    samples = set()
    with open(samples_csv_path, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "sample" not in reader.fieldnames:
            raise ValueError(f"Column 'sample' missing in {samples_csv_path}. Found: {reader.fieldnames}")
        for row in reader:
            val = (row.get("sample") or "").strip()
            if val:
                samples.add(val)

    known = sorted(samples, key=len, reverse=True)
    print(f"[INFO] Loaded {len(known)} sample IDs from {samples_csv_path}")
    return known


def extract_sample_id(header, known_samples):
    """
    Extracts a sample ID from a FASTA header.
    Matches exact token or prefix followed by a common delimiter ('_', '|', ':', '/', ';', ',').
    """
    token = header.strip().split()[0] if header.strip() else ""
    delimiters = ("_", "|", ":", "/", ";", ",")

    for sample in known_samples:
        if token == sample:
            return sample
        for delim in delimiters:
            if token.startswith(sample + delim):
                return sample

    if known_samples:
        print(f"[WARN] Header token '{token}' did not match any known sample; keeping raw token.", file=sys.stderr)
    return token


def parse_fasta_seqs(input_fa, known_samples):
    """
    Parses a FASTA file into a dict mapping sample_id -> uppercase aligned sequence.
    """
    seqs = {}
    curr_sample = None
    curr_seq = []

    def _commit(sample_id, parts):
        if sample_id is None or not parts:
            return
        seq_str = "".join(parts).upper()
        if sample_id in seqs and seqs[sample_id] != seq_str:
            print(f"[WARN] Duplicate sample '{sample_id}' in {input_fa}; overwriting previous sequence.", file=sys.stderr)
        seqs[sample_id] = seq_str

    with open(input_fa, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _commit(curr_sample, curr_seq)
                curr_sample = extract_sample_id(line[1:], known_samples)
                curr_seq = []
            else:
                curr_seq.append(line)

        _commit(curr_sample, curr_seq)

    return seqs


def analyze_gene_alignment(gene_name, seqs, standard_samples):
    """
    Analyzes a single alignment for length, variable sites, parsimony-informative sites, and gaps.
    """
    max_len = max((len(s) for s in seqs.values()), default=0)
    ordered_seqs = [seqs.get(s, "").ljust(max_len, "-")[:max_len] for s in standard_samples]

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

        informative_chars = [c for c in col_chars if c not in NON_INFORMATIVE]
        unique_states = set(informative_chars)

        if len(unique_states) > 1:
            variable_sites += 1
            multiton_states = sum(1 for state in unique_states if informative_chars.count(state) >= 2)
            if multiton_states >= 2:
                parsimony_informative += 1

    denom = num_taxa * aln_length
    return {
        "gene": gene_name,
        "taxa_count": num_taxa,
        "alignment_length": aln_length,
        "variable_sites": variable_sites,
        "parsimony_informative_sites": parsimony_informative,
        "gap_columns": gap_columns,
        "total_gaps": total_gaps,
        "fraction_gaps": round(total_gaps / denom, 4) if denom > 0 else 0,
    }


def write_normalized_fasta(seqs, output_fa, standard_samples):
    """
    Writes a normalized FASTA with all standard_samples padded to max_len using '-'.
    """
    max_len = max((len(s) for s in seqs.values()), default=0)
    if max_len == 0:
        raise ValueError(f"Empty alignment for {output_fa} (max_len=0).")

    with open(output_fa, "w", encoding="utf-8") as out:
        for sample in standard_samples:
            padded = seqs.get(sample, "").ljust(max_len, "-")[:max_len]
            out.write(f">{sample}\n{padded}\n")


def extract_gene_name(filename):
    base_name = os.path.basename(filename)
    for suffix in ("_trimmed.fa", "_trimmed.fasta", ".fasta", ".fa", ".fna"):
        if base_name.endswith(suffix):
            return base_name[: -len(suffix)]
    return base_name


def parse_args():
    parser = argparse.ArgumentParser(description="Normalize gene alignments, audit variation, and concatenate via PhyKIT.")
    parser.add_argument("--samples-table", default=None, help="Path to samples CSV containing a 'sample' column.")
    parser.add_argument("alignments", nargs="+", help="Input alignment FASTA files.")
    return parser.parse_args()


def main():
    args = parse_args()
    known_samples = load_known_samples(args.samples_table)

    parsed_alignments = []
    all_samples = set(known_samples)

    for fa in args.alignments:
        if not os.path.isfile(fa):
            print(f"[ERROR] Input file not found: {fa}", file=sys.stderr)
            sys.exit(1)
        seqs = parse_fasta_seqs(fa, known_samples)
        if not seqs:
            print(f"[ERROR] No sequences found in {fa}", file=sys.stderr)
            sys.exit(1)
        parsed_alignments.append((fa, seqs))
        all_samples.update(seqs.keys())

    if not all_samples:
        print("[ERROR] No samples found in table or alignments.", file=sys.stderr)
        sys.exit(1)

    standard_samples = sorted(all_samples)
    print(f"[INFO] Fixed sample order ({len(standard_samples)} taxa): {standard_samples}")

    cleaned_dir = "cleaned_for_phykit"
    os.makedirs(cleaned_dir, exist_ok=True)

    gene_stats_list = []
    cleaned_list_file = "alignment_list_cleaned.txt"

    with open(cleaned_list_file, "w", encoding="utf-8") as list_out:
        for fa, seqs in parsed_alignments:
            base_name = os.path.basename(fa)
            gene_name = extract_gene_name(base_name)
            out_fa = os.path.join(cleaned_dir, base_name)

            write_normalized_fasta(seqs, out_fa, standard_samples)
            gene_stats_list.append(analyze_gene_alignment(gene_name, seqs, standard_samples))
            list_out.write(os.path.abspath(out_fa) + "\n")

    report_file = "gene_variation_report.csv"
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
    with open(report_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(gene_stats_list)
    print(f"[SUCCESS] Generated {report_file}")

    cmd = ["phykit", "create_concat", "-a", cleaned_list_file, "-p", "concatenated_output"]
    print(f"[INFO] Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] PhyKIT failed: {e}", file=sys.stderr)
        sys.exit(1)

    if os.path.exists("concatenated_output.fa"):
        os.replace("concatenated_output.fa", "concatenated_alignment.fasta")
    if os.path.exists("concatenated_output.partition"):
        os.replace("concatenated_output.partition", "partitions.txt")


if __name__ == "__main__":
    main()
