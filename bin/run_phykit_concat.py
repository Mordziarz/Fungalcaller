#!/usr/bin/env python3
import os
import sys
import subprocess
import csv

def analyze_gene_alignment(gene_name, seqs, standard_samples):
    """
    Analyzes a single gene alignment for length, mutations (variability), and gaps (indels).
    """
    ordered_seqs = [seqs.get(s, "") for s in standard_samples]
    
    max_len = max(len(s) for s in ordered_seqs) if ordered_seqs else 0
    ordered_seqs = [s.ljust(max_len, '-') for s in ordered_seqs]
    
    num_taxa = len(standard_samples)
    aln_length = max_len
    
    variable_sites = 0
    parsimony_informative = 0
    gap_columns = 0
    total_gaps = 0

    for col_idx in range(aln_length):
        col_chars = [seq[col_idx] for seq in ordered_seqs]
        
        if '-' in col_chars:
            gap_columns += 1
            total_gaps += col_chars.count('-')
            
        non_gap_chars = [c for c in col_chars if c != '-']
        
        if len(set(non_gap_chars)) > 1:
            variable_sites += 1
            char_counts = {c: non_gap_chars.count(c) for c in set(non_gap_chars)}
            multiton_count = sum(1 for c, count in char_counts.items() if count >= 2)
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
        "fraction_gaps": round(total_gaps / (num_taxa * aln_length), 4) if (num_taxa * aln_length) > 0 else 0
    }

def process_fasta_for_concatenation(input_fa, output_fa, standard_samples):
    seqs = {}
    with open(input_fa, "r") as f:
        curr_sample = None
        curr_seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_sample and curr_seq:
                    seqs[curr_sample] = "".join(curr_seq)
                header = line[1:]
                matched = next((s for s in standard_samples if s in header), header.split("_")[0])
                curr_sample = matched
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_sample and curr_seq:
            seqs[curr_sample] = "".join(curr_seq)

    with open(output_fa, "w") as out:
        for sample in standard_samples:
            if sample in seqs:
                out.write(f">{sample}\n{seqs[sample]}\n")
            else:
                out.write(f">{sample}\n" + "-" * len(next(iter(seqs.values()))) + "\n")

    return seqs

def main():
    input_files = sys.argv[1:]
    if not input_files:
        print("[ERROR] No input files provided for PhyKIT!", file=sys.stderr)
        sys.exit(1)

    all_samples = set()
    for fa in input_files:
        with open(fa, "r") as f:
            for line in f:
                if line.startswith(">"):
                    h = line[1:].strip()
                    s_name = h.split("_")[0]
                    all_samples.add(s_name)
    
    standard_samples = sorted(list(all_samples))
    print(f"[INFO] Established fixed sample sorting order: {standard_samples}")

    cleaned_dir = "cleaned_for_phykit"
    os.makedirs(cleaned_dir, exist_ok=True)

    gene_stats_list = []
    cleaned_list_file = "alignment_list_cleaned.txt"

    with open(cleaned_list_file, "w") as list_out:
        for fa in input_files:
            base_name = os.path.basename(fa)
            gene_name = base_name.replace("_trimmed.fa", "").replace(".fa", "")
            out_fa = os.path.join(cleaned_dir, base_name)
            
            seqs = process_fasta_for_concatenation(fa, out_fa, standard_samples)
            
            stats = analyze_gene_alignment(gene_name, seqs, standard_samples)
            gene_stats_list.append(stats)
            
            list_out.write(os.path.abspath(out_fa) + "\n")

    report_file = "gene_variation_report.csv"
    with open(report_file, "w", newline="") as csvfile:
        fieldnames = [
            "gene", "taxa_count", "alignment_length", 
            "variable_sites", "parsimony_informative_sites", 
            "gap_columns", "total_gaps", "fraction_gaps"
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
        os.rename("concatenated_output.fa", "concatenated_alignment.fasta")
    if os.path.exists("concatenated_output.partition"):
        os.rename("concatenated_output.partition", "partitions.txt")
    if os.path.exists(report_file):
        pass

if __name__ == "__main__":
    main()
