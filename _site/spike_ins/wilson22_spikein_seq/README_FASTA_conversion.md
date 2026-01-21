# Table S1 to FASTA Conversion - README

## Overview

This directory contains scripts to convert the Wilson et al. 2022 Table S1 (spike-in control sequences) from Excel format to FASTA format.

## Files Included

1. **table_s1_to_fasta.R** - Full R script with detailed output
2. **table_s1_to_fasta_simple.R** - Minimal R code (10 lines)
3. **table_s1_to_fasta.py** - Python script version
4. **spike_in_controls.fa** - Pre-generated FASTA file (52 sequences)

## Input File Structure

**Table S1** contains 52 spike-in control sequences with:
- **name**: Sequence identifier (e.g., "80b_1C_35G-1")
- **seq**: DNA sequence
- **methylated**: Methylation status (TRUE/FALSE)
- **len**: Fragment length (80, 160, or 320 bp)
- **gc**: G+C content (35%, 50%, or 65%)
- **cpg_frac**: CpG fraction (0.0125, 0.025, or 0.05)

### Sequence Distribution

- **80 bp fragments**: 18 sequences (9 methylated, 9 unmethylated)
- **160 bp fragments**: 16 sequences (8 methylated, 8 unmethylated)
- **320 bp fragments**: 18 sequences (9 methylated, 9 unmethylated)
- **Total**: 52 sequences (26 methylated, 26 unmethylated)

## Usage

### Option 1: Use Pre-generated FASTA File

The file `spike_in_controls.fa` is already generated and ready to use.

```bash
# Use directly in your analysis
bowtie2-build spike_in_controls.fa spike_in_controls
```

### Option 2: Run R Script (Full Version)

```bash
# Make sure TableS1_spike_in_sequences.xlsx is in the same directory
Rscript table_s1_to_fasta.R
```

This will:
- Read the Excel file
- Display summary statistics
- Create `spike_in_controls.fa`
- Show preview of output

### Option 3: Run R Script (Simple Version)

```R
# In R console or RStudio
source("table_s1_to_fasta_simple.R")
```

Or from command line:
```bash
Rscript table_s1_to_fasta_simple.R
```

### Option 4: Run Python Script

```bash
# Make sure pandas and openpyxl are installed
pip install pandas openpyxl

# Run the script
python table_s1_to_fasta.py
```

## FASTA Format

The output FASTA file has standard format:

```
>80b_1C_35G-1
TAGGATATAGGTTGTCCCCTAGTAGGAGATAAACTTTGATTAACATCCAATTGATCGTTAGTGTCCTTCAAAATTATGCT
>80b_1C_35G-2
TGTCTAAATTAAAGTTGTGATCTTTGACTTAGCAACGTCTCACCCTATAGCCTACCAGACAAGAATTATGAAGAACATAT
...
```

Each sequence entry has:
- Header line starting with `>` followed by sequence name
- Sequence line with DNA bases (A, T, G, C)

## Using the FASTA File

### Create Bowtie2 Index

```bash
# Build index for alignment
bowtie2-build spike_in_controls.fa spike_in_controls

# This creates index files:
# - spike_in_controls.1.bt2
# - spike_in_controls.2.bt2
# - spike_in_controls.3.bt2
# - spike_in_controls.4.bt2
# - spike_in_controls.rev.1.bt2
# - spike_in_controls.rev.2.bt2
```

### Align Reads to Spike-ins

```bash
# Align paired-end reads
bowtie2 -x spike_in_controls \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -S spike_in_alignment.sam

# Count aligned reads
samtools view -c -F 4 spike_in_alignment.sam
```

### Extract Spike-in Reads

```bash
# Get reads that mapped to spike-ins
samtools view -b -F 4 spike_in_alignment.sam > spike_in_mapped.bam

# Get reads that did NOT map to spike-ins (for alignment to human genome)
samtools view -b -f 4 spike_in_alignment.sam > unmapped.bam
samtools fastq -1 unmapped_R1.fq -2 unmapped_R2.fq unmapped.bam
```

## Sequence Naming Convention

Sequence names follow the pattern: `{length}b_{CpGs}C_{GC}G-{replicate}`

Examples:
- `80b_1C_35G-1`: 80 bp, 1 CpG per 80 bp (0.0125 CpG fraction), 35% G+C, replicate 1
- `160b_4C_50G-2`: 160 bp, 4 CpG per 40 bp (0.025 CpG fraction), 50% G+C, replicate 2
- `320b_16C_65G-1`: 320 bp, 16 CpG per 20 bp (0.05 CpG fraction), 65% G+C, replicate 1

Methylated vs unmethylated versions have different sequences but same identifiers.

## Verification

### Check File Contents

```bash
# Count sequences
grep -c "^>" spike_in_controls.fa
# Should output: 52

# Check sequence lengths
grep -v "^>" spike_in_controls.fa | awk '{print length($0)}' | sort | uniq -c
# Should show: 18 sequences of 80bp, 16 of 160bp, 18 of 320bp
```

### R Verification

```R
library(Biostrings)
seqs <- readDNAStringSet("spike_in_controls.fa")
length(seqs)  # Should be 52
table(width(seqs))  # Shows length distribution
```

### Python Verification

```python
from Bio import SeqIO

seqs = list(SeqIO.parse("spike_in_controls.fa", "fasta"))
print(f"Total sequences: {len(seqs)}")
lengths = [len(seq) for seq in seqs]
print(f"Length distribution: {set(lengths)}")
```

## Troubleshooting

### Missing readxl package (R)

```R
install.packages("readxl")
```

### Missing openpyxl package (Python)

```bash
pip install openpyxl pandas
```

### Excel file not found

Make sure `TableS1_spike_in_sequences.xlsx` is in the same directory as the script.

## References

Wilson SL, Shen SY, Harmon L, et al. Sensitive and reproducible cell-free methylome quantification with synthetic spike-in controls. *Cell Rep Methods*. 2022;2(9):100294. doi:10.1016/j.crmeth.2022.100294

---

*Generated: January 2026*
