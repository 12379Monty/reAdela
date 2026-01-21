# Data Download and Organization Guide: Wilson et al. 2022 Spike-in Paper

## Paper Citation
Wilson SL, Shen SY, Harmon L, Burgener JM, Triche T Jr, Bratman SV, De Carvalho DD, Hoffman MM. Sensitive and reproducible cell-free methylome quantification with synthetic spike-in controls. *Cell Rep Methods*. 2022 Sep 9;2(9):100294. doi: 10.1016/j.crmeth.2022.100294. PMID: 36160046.

---

## Repository and Code Locations

### Main Code Repository
- **GitHub**: https://github.com/hoffmangroup/2020spikein
- **Zenodo** (archived version): https://doi.org/10.5281/zenodo.4683791
- **GitHub Pages** (documentation): https://hoffmangroup.github.io/2020spikein/

### R Package
- **Bioconductor spiky package**: https://bioconductor.org/packages/spiky
- **GitHub**: https://github.com/trichelab/spiky

---

## Data Availability Summary

| Data Type | Location | Access | Accession |
|-----------|----------|--------|-----------|
| HCT116 cell line cfMeDIP-seq | GEO | Public | GSE166259 |
| AML patient raw data | EGA | Controlled | EGAS00001005069 |
| AML patient processed data | GEO | Public | GSE166259 |
| Code and scripts | GitHub/Zenodo | Public | github.com/hoffmangroup/2020spikein |

---

## Part 1: Setting Up Your Working Directory

### Recommended Directory Structure

```
wilson2022_reproduction/
├── code/
│   └── 2020spikein/          # Cloned GitHub repository
├── data/
│   ├── raw/
│   │   ├── HCT116/           # HCT116 cell line data from GEO
│   │   └── AML/              # AML patient data (if accessible)
│   ├── processed/            # Intermediate processed files
│   └── spike_in_sequences/   # Spike-in control sequences
├── results/
│   ├── figures/              # Reproduced figures
│   └── tables/               # Reproduced tables
└── reference/
    ├── genome/               # Human reference genome (hg38)
    └── annotations/          # ENCODE blacklist, RepeatMasker, etc.
```

### Create Directory Structure

```bash
mkdir -p wilson2022_reproduction/{code,data/{raw/{HCT116,AML},processed,spike_in_sequences},results/{figures,tables},reference/{genome,annotations}}
cd wilson2022_reproduction
```

---

## Part 2: Downloading the Code

### Option 1: Clone from GitHub (if network permits)

```bash
cd code
git clone https://github.com/hoffmangroup/2020spikein.git
```

### Option 2: Download from Zenodo

1. Visit: https://doi.org/10.5281/zenodo.4683791
2. Download the archive: `2020spikein-2020spikein-v.1.0.0.zip`
3. Extract to `code/2020spikein/`

```bash
# If downloaded to Downloads folder
unzip ~/Downloads/2020spikein-2020spikein-v.1.0.0.zip -d code/
mv code/swils6-2020spikein-* code/2020spikein
```

### Repository Contents

The repository contains:
- `Preprocessing/` - Scripts for processing raw FASTQ files
- `Analysis/` - Scripts for analysis and figure generation
- `figures/` - Figure outputs
- `README.md` - Basic documentation
- `index.Rmd` and `index.html` - Detailed documentation

---

## Part 3: Downloading Sequencing Data from GEO

### Access GEO Dataset GSE166259

**Direct link**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166259

This dataset contains:
1. HCT116 cell line cfMeDIP-seq data (technical replicates)
2. AML patient processed data (molar amount estimates)

### Download Methods

#### Method 1: Using SRA Toolkit (Recommended for FASTQ files)

```bash
# Install SRA Toolkit if not already installed
# https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

# Download HCT116 samples
# First, get the SRA run accessions from GEO page

# Example for downloading one sample:
prefetch SRR1234567  # Replace with actual SRR number
fastq-dump --split-files --gzip SRR1234567 -O data/raw/HCT116/

# For parallel download of multiple samples:
cat > data/raw/sra_accessions.txt << EOF
SRR1234567
SRR1234568
# Add all SRR accessions here
EOF

# Batch download
parallel -j 4 'prefetch {} && fastq-dump --split-files --gzip {} -O data/raw/HCT116/' :::: data/raw/sra_accessions.txt
```

#### Method 2: Using GEO FTP (for supplementary files)

```bash
# Access GEO FTP site
wget -r -np -nd -P data/raw/HCT116/ \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/suppl/

# Alternatively, download specific files from GEO web interface
# and save to data/raw/HCT116/
```

#### Method 3: Using NCBI Datasets Command-line Tool

```bash
# Install datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

datasets download genome accession GSE166259 \
  --include gff3,rna,cds,protein,genome,seq-report \
  --filename data/raw/GSE166259.zip

unzip data/raw/GSE166259.zip -d data/raw/HCT116/
```

### Expected HCT116 Data Files

Based on the paper, you should have:
- **Input samples**: 10% of DNA before cfMeDIP (2 technical replicates)
- **Output samples**: After cfMeDIP-seq (2 technical replicates)
- **Spike-in amounts**: 0.01 ng, 0.05 ng, 0.1 ng spike-in DNA + 10 ng HCT116
- **EPIC array data**: HCT116 methylation array data for validation

---

## Part 4: Downloading AML Patient Data

### Controlled Access Data (EGA)

**Study**: EGAS00001005069  
**Link**: https://ega-archive.org/studies/EGAS00001005069

#### Access Requirements

1. **Create EGA Account**: https://ega-archive.org/register
2. **Request Access**:
   - Contact: University Health Network Genomics Data Access Committee
   - Email: dac@uhn.ca
   - Policy: https://doi.org/10.5281/zenodo.4568265

3. **Application Process**:
   - Submit data access request form
   - Provide research justification
   - Sign data use agreement
   - Typical approval time: 2-4 weeks

#### After Approval: Download from EGA

```bash
# Install pyEGA3 client
pip install pyega3 --break-system-packages

# Login to EGA
pyega3 -cf path/to/credentials_file.json

# List datasets
pyega3 -cf path/to/credentials_file.json datasets

# Download dataset
pyega3 -cf path/to/credentials_file.json fetch EGAD00001007XXX --saveto data/raw/AML/

# Decrypt files (if necessary)
# Follow EGA-specific decryption instructions
```

### Alternative: Use Processed AML Data from GEO

The processed AML data (molar amount estimates) is publicly available on GEO:GSE166259

```bash
# Download processed data
wget -P data/processed/ \
  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE166259&format=file

# Extract
tar -xzf data/processed/GSE166259_RAW.tar -C data/processed/
```

---

## Part 5: Spike-in Control Sequences

### Download Spike-in Sequences

The 54 spike-in control sequences are documented in **Table S1** of the paper.

#### Option 1: Extract from Supplementary Materials

1. Visit paper on PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC9499995/
2. Download **Table S1** (mmc1.xlsx)
3. Save to `data/spike_in_sequences/`

```bash
# If downloaded
mv ~/Downloads/mmc1.xlsx data/spike_in_sequences/TableS1_spike_in_sequences.xlsx
```

#### Option 2: From GitHub Repository

The spike-in sequences should be in the repository or can be extracted from the scripts.

```bash
# Check if sequences are in repository
grep -r "spike.*in.*sequence" code/2020spikein/
```

#### Option 3: Manual Entry from Paper

The paper describes the spike-in design:
- **Fragment lengths**: 80 bp, 160 bp, 320 bp
- **G+C content**: 35%, 50%, 65%
- **CpG fraction**: 1/80 bp, 1/40 bp, 1/20 bp
- **Total**: 52 fragments (2 failed synthesis)
- **Methylation status**: Methylated and unmethylated versions

### Create Spike-in Reference FASTA

```bash
# Create a FASTA file with spike-in sequences
# This will be needed for alignment

cat > data/spike_in_sequences/spikein_controls.fa << 'EOF'
>spikein_80bp_35GC_1per80_meth_1
[sequence from Table S1]
>spikein_80bp_35GC_1per80_unmeth_1
[sequence from Table S1]
# ... add all 52 sequences
EOF
```

---

## Part 6: Reference Genome and Annotations

### Human Reference Genome (GRCh38/hg38)

```bash
cd reference/genome

# Download from UCSC
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Or download from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38.fa

# Index for Bowtie2
bowtie2-build hg38.fa hg38
```

### ENCODE Blacklist

```bash
cd reference/annotations

# Download ENCODE blacklist for hg38
wget -O hg38_blacklist.bed.gz \
  https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38_blacklist.bed.gz
```

### RepeatMasker

```bash
# Download RepeatMasker annotations
wget -O hg38_rmsk.txt.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
gunzip hg38_rmsk.txt.gz

# Process into BED format
awk 'BEGIN{OFS="\t"}{print $6,$7,$8,$11,$12,$13}' hg38_rmsk.txt > hg38_rmsk.bed
```

### UCSC Simple Repeats

```bash
# Download simple repeats
wget -O hg38_simpleRepeat.txt.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
gunzip hg38_simpleRepeat.txt.gz

# Convert to BED
awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' hg38_simpleRepeat.txt > hg38_simpleRepeat.bed
```

### Umap Mappability Scores

```bash
# Download Umap k100 mappability
wget https://bismap.hoffmanlab.org/data/hg38/k100.umap.bedgraph.gz
gunzip k100.umap.bedgraph.gz
```

### GENCODE Annotations

```bash
# Download GENCODE v33 (used in paper)
wget -O gencode.v33.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
gunzip gencode.v33.annotation.gtf.gz
```

---

## Part 7: Installing Required Software

### R and Bioconductor Packages

```R
# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install spiky package
BiocManager::install("spiky")

# Install other required packages
required_packages <- c(
  "sesame",          # For EPIC array data
  "compute.es",      # For effect size calculations
  "GenomicRanges",   # For genomic intervals
  "rtracklayer",     # For importing/exporting genomic data
  "ggplot2",         # For plotting
  "dplyr",           # For data manipulation
  "tidyr"            # For data tidying
)

BiocManager::install(required_packages)
```

### Command-line Tools

```bash
# fastp (v0.11.5 or later)
# https://github.com/OpenGene/fastp
conda install -c bioconda fastp

# Bowtie2 (v2.3.5 or later)
# http://bowtie-bio.sourceforge.net/bowtie2/
conda install -c bioconda bowtie2

# SAMtools (v0.10.2 or later)
# https://github.com/samtools/samtools
conda install -c bioconda samtools

# UMI-tools (v1.0.0 or later)
# https://github.com/CGATOxford/UMI-tools
pip install umi_tools --break-system-packages

# bedtools (v2.29.2 or later)
# https://github.com/arq5x/bedtools2
conda install -c bioconda bedtools
```

---

## Part 8: Data Organization Checklist

### Before Running Scripts

Ensure you have the following files organized:

#### Raw Sequencing Data
- [ ] HCT116 FASTQ files (input and output, multiple spike-in amounts)
- [ ] AML patient FASTQ files (3 labs × 5 samples) OR processed data

#### Spike-in Reference
- [ ] Spike-in sequences FASTA file
- [ ] Table S1 with spike-in properties (fragment length, GC, CpG fraction)
- [ ] Bowtie2 index for spike-in sequences

#### Human Reference
- [ ] hg38 reference genome FASTA
- [ ] Bowtie2 index for hg38
- [ ] ENCODE blacklist BED file
- [ ] RepeatMasker BED file
- [ ] Simple repeats BED file
- [ ] Umap k100 mappability scores
- [ ] GENCODE v33 annotations

#### Array Data (if reproducing array comparisons)
- [ ] HCT116 EPIC array data (IDAT files or processed)

---

## Part 9: Running the Analysis Pipeline

### Step 1: Preprocessing

The `Preprocessing/` folder contains numbered scripts:

```bash
cd code/2020spikein/Preprocessing

# 1. Trim adapters and QC
# 2. Align to spike-in reference
# 3. Align unaligned reads to hg38
# 4. Remove duplicates using UMIs
# 5. Generate read count matrices
```

### Step 2: Analysis

The `Analysis/` folder contains:

```bash
cd code/2020spikein/Analysis

# 1. Train GLM model on spike-in data
# 2. Predict molar amounts for genomic windows
# 3. Compare with EPIC array data
# 4. Assess batch effects
# 5. Generate figures
```

### Step 3: Run spiky Package

```R
library(spiky)

# Load your processed data
# Follow vignettes at:
# https://bioconductor.org/packages/release/bioc/vignettes/spiky/inst/doc/spiky.html
```

---

## Part 10: Expected Outputs

### Figures to Reproduce

From the paper:
- **Figure 1**: Experimental design
- **Figure 2**: Assessing biases (fragment length, G+C, CpG)
- **Figure 3**: 2D histograms (molar amount vs. SD and mappability)
- **Figure 4**: Correlation with EPIC array M-values
- **Figure 5**: Model accuracy (held-out spike-ins)
- **Figure 6**: PCA and batch effects

### Tables to Reproduce

- **Table 1**: Range of reads for alternative fragments
- **Table 2**: High molar amount genomic windows
- **Table S2-S6**: Supplementary data tables

---

## Troubleshooting

### Issue: Cannot Access GitHub Repository

**Solution**: Download from Zenodo (https://doi.org/10.5281/zenodo.4683791)

### Issue: EGA Access Denied

**Solution**: 
1. Use processed data from GEO (GSE166259) instead
2. Focus on HCT116 data which is publicly available
3. Contact dac@uhn.ca with research justification

### Issue: Missing Spike-in Sequences

**Solution**: 
1. Download Table S1 from PMC article
2. Extract sequences from supplementary materials
3. Contact paper authors if sequences not available

### Issue: Software Version Incompatibilities

**Solution**:
- Check exact versions in paper's Methods section
- Use conda environments to manage versions
- Refer to original protocol papers for alternative tools

---

## Additional Resources

### Paper Links
- **PubMed**: https://pubmed.ncbi.nlm.nih.gov/36160046/
- **PMC**: https://pmc.ncbi.nlm.nih.gov/articles/PMC9499995/
- **Cell Reports Methods**: https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(22)00176-X

### Protocol Papers
- **cfMeDIP-seq Protocol** (Shen et al. 2019): https://www.nature.com/articles/s41596-019-0202-2
- **Original cfMeDIP-seq** (Shen et al. 2018): https://www.nature.com/articles/s41586-018-0703-0

### Data Repositories
- **GEO**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166259
- **EGA**: https://ega-archive.org/studies/EGAS00001005069
- **Zenodo**: https://doi.org/10.5281/zenodo.4683791

### Contact
- **Lead Author**: Michael M. Hoffman (michael.hoffman@utoronto.ca)
- **EGA Data Access**: dac@uhn.ca

---

## Quick Start Summary

For a minimal reproduction focused on HCT116 data:

```bash
# 1. Create directory structure
mkdir -p wilson2022/{code,data,results,reference}
cd wilson2022

# 2. Download code
wget https://zenodo.org/record/4683791/files/2020spikein.zip
unzip 2020spikein.zip -d code/

# 3. Download HCT116 data from GEO
# Use SRA toolkit or web interface for GSE166259

# 4. Download spike-in sequences
# From paper supplementary materials (Table S1)

# 5. Download reference genome and annotations
# Use commands in Part 6

# 6. Install required software
# R packages and command-line tools from Part 7

# 7. Run preprocessing scripts
cd code/2020spikein/Preprocessing
# Follow numbered scripts in order

# 8. Run analysis scripts
cd ../Analysis
# Follow numbered scripts in order
```

---

*Last Updated: January 2026*
