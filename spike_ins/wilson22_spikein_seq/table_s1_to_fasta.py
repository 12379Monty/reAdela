#!/usr/bin/env python3

"""
Python script to convert Table S1 (spike-in sequences) to FASTA format
Wilson et al. 2022 - Cell Reports Methods
"""

import pandas as pd
import sys

def write_fasta(df, output_file):
    """
    Write sequences to FASTA format
    
    Args:
        df: DataFrame with 'name' and 'seq' columns
        output_file: Path to output FASTA file
    """
    with open(output_file, 'w') as f:
        for idx, row in df.iterrows():
            # Write FASTA header (starts with >)
            f.write(f">{row['name']}\n")
            # Write sequence
            f.write(f"{row['seq']}\n")
    
    print(f"FASTA file written to: {output_file}")
    print(f"Number of sequences: {len(df)}")

def main():
    # Input and output file paths
    input_file = "TableS1_spike_in_sequences.xlsx"
    output_file = "spike_in_controls.fa"
    
    print(f"Reading Excel file: {input_file}")
    
    try:
        # Read the Excel file - Table S1 sheet contains the actual data
        df = pd.read_excel(input_file, sheet_name="Table S1")
    except FileNotFoundError:
        print(f"ERROR: File '{input_file}' not found!")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR reading file: {e}")
        sys.exit(1)
    
    # Display information about the data
    print(f"\nData summary:")
    print(f"  Total sequences: {len(df)}")
    print(f"  Columns: {', '.join(df.columns)}")
    
    # Show length distribution
    print(f"\nSequence length distribution:")
    print(df['len'].value_counts().sort_index())
    
    # Show methylation status
    print(f"\nMethylation status:")
    print(df['methylated'].value_counts())
    
    # Verify required columns exist
    if 'name' not in df.columns or 'seq' not in df.columns:
        print("ERROR: Required columns 'name' and/or 'seq' not found!")
        sys.exit(1)
    
    # Check for missing data
    if df['name'].isna().any() or df['seq'].isna().any():
        print("ERROR: Missing data in name or sequence columns!")
        sys.exit(1)
    
    # Write to FASTA file
    print(f"\nWriting FASTA file...")
    write_fasta(df, output_file)
    
    # Show preview
    print(f"\nFirst 6 lines of {output_file}:")
    with open(output_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= 6:
                break
            print(line.rstrip())
    
    print("\nDone!")

if __name__ == "__main__":
    main()
