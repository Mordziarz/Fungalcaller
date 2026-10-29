#!/usr/bin/env python3
import os
import glob
import pandas as pd
import sys

def summarize_jellyfish(input_dir, output_csv):
    histo_files = glob.glob(os.path.join(input_dir, "*_histo.txt"))
    
    if not histo_files:
        print(f"[WARNING] No *_histo.txt files found in directory: {input_dir}")
        return

    combined_df = None

    for file_path in histo_files:
        base_name = os.path.basename(file_path)
        sample_name = base_name.replace("_histo.txt", "")
        
        df = pd.read_csv(file_path, sep=r'\s+', header=None, names=['Coverage', sample_name])
        
        if combined_df is None:
            combined_df = df
        else:
            combined_df = pd.merge(combined_df, df, on='Coverage', how='outer')

    if combined_df is not None:
        combined_df['Coverage'] = pd.to_numeric(combined_df['Coverage'], errors='coerce')
        combined_df = combined_df.dropna(subset=['Coverage']).sort_values(by='Coverage')
        
        combined_df = combined_df.fillna(0)
        
        for col in combined_df.columns:
            combined_df[col] = combined_df[col].astype(int)

        combined_df.to_csv(output_csv, index=False)
        print(f"[SUCCESS] Generated consolidated histogram summary table: {output_csv}")

if __name__ == "__main__":
    input_directory = sys.argv[1] if len(sys.argv) > 1 else "."
    summarize_jellyfish(input_directory, "jellyfish_summary_table.csv")
