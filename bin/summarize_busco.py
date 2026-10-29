#!/usr/bin/env python3
import os
import glob
import csv
import re
import sys

def parse_busco_summaries(results_dir, output_file="busco_summary_table.csv"):
    summary_files = glob.glob(os.path.join(results_dir, "**/short_summary.*.txt"), recursive=True)
    
    table_data = []
    
    for file_path in summary_files:
        parts = file_path.split(os.sep)
        raw_sample_name = parts[-2] if len(parts) >= 2 else "unknown"
        sample_name = raw_sample_name.replace("_busco_out", "")
        
        c_pct, s_pct, d_pct, f_pct, m_pct, total = "", "", "", "", "", ""
        
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if "C:" in line and "%" in line and "n:" in line:
                    match = re.search(r'C:([\d\.]+)%\[S:([\d\.]+)%,D:([\d\.]+)%\]\s*,F:([\d\.]+)%,\s*M:([\d\.]+)%,\s*n:(\d+)', line)
                    if match:
                        c_pct, s_pct, d_pct, f_pct, m_pct, total = match.groups()
                        
        if c_pct:
            table_data.append({
                'Sample': sample_name,
                'Complete (%)': c_pct,
                'Single-copy (%)': s_pct,
                'Duplicated (%)': d_pct,
                'Fragmented (%)': f_pct,
                'Missing (%)': m_pct,
                'Total BUSCOs': total
            })
            
    if table_data:
        table_data = sorted(table_data, key=lambda x: x['Sample'])
        keys = table_data[0].keys()
        with open(output_file, 'w', newline='', encoding='utf-8') as output_csv:
            writer = csv.DictWriter(output_csv, fieldnames=keys)
            writer.writeheader()
            writer.writerows(table_data)
        print(f"[SUCCESS] Generated summary table: {output_file}")
    else:
        print(f"[WARNING] No BUSCO summary files found in: {results_dir}")

if __name__ == "__main__":
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "./results/busco"
    parse_busco_summaries(target_dir)
