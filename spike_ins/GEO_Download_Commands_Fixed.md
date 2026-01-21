# Corrected Commands for Downloading GEO Data - GSE166259

## The Problem
The original command had an unquoted URL with special characters, which causes shell parsing issues.

## Corrected Methods for Downloading from GEO

### Method 1: Direct FTP Download (Recommended)

```bash
# Download all supplementary files from the series
# GSE166259 FTP structure: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/

# Download supplementary files
wget -r -np -nd -P data/raw/GSE166259/ \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/suppl/

# Or download the entire series directory
wget -r -np -P data/raw/GSE166259/ \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/
```

### Method 2: Using curl

```bash
# Download supplementary files with curl
curl -o data/raw/GSE166259_suppl.tar \
  'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/suppl/GSE166259_RAW.tar'

# Extract the tarball
tar -xvf data/raw/GSE166259_suppl.tar -C data/raw/GSE166259/
```

### Method 3: Using wget with quoted URL

```bash
# If you need to use the HTTP interface, quote the URL properly
wget -O data/raw/GSE166259_RAW.tar \
  'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE166259&format=file'

# Or using double quotes with proper escaping
wget -O data/raw/GSE166259_RAW.tar \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE166259&format=file"
```

### Method 4: Using R/Bioconductor GEOquery (Most Reliable)

This is often the easiest and most reliable method:

```R
# Install GEOquery if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

# Set timeout for large downloads
options(timeout = 600)

# Download supplementary files
getGEOSuppFiles("GSE166259", baseDir = "data/raw/")

# This will create data/raw/GSE166259/ directory with all supplementary files

# Alternatively, get the full GEO object
gse <- getGEO("GSE166259", GSEMatrix = TRUE, destdir = "data/raw/")
```

### Method 5: Using Python GEOparse

```python
import GEOparse

# Download the series
gse = GEOparse.get_GEO(geo="GSE166259", destdir="./data/raw/")

# Download supplementary files
gse.download_supplementary_files(directory="./data/raw/GSE166259/", 
                                  download_sra=False)
```

### Method 6: Manual Download from Web Browser

1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166259
2. Scroll to bottom of page
3. Click on "Download family" section
4. Click link for supplementary files or series matrix files
5. Save to `data/raw/GSE166259/`

## Understanding GEO File Structure

### Files You'll Find:

1. **SOFT files** (`GSE166259_family.soft.gz`):
   - Contains all metadata
   - Text format with structured information

2. **Series Matrix** (`GSE166259_series_matrix.txt.gz`):
   - Processed data in matrix format
   - Sample metadata

3. **Supplementary Files** (in `suppl/` directory):
   - Raw data files (FASTQ, BAM, etc.)
   - Processed data files (count matrices, normalized values)
   - Sample-specific files (one per GSM accession)

4. **RAW tarball** (`GSE166259_RAW.tar`):
   - Contains all supplementary files
   - Convenient single download

## Checking What's Available

### Using curl to list FTP directory:

```bash
# List supplementary files
curl -l ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/suppl/

# Or with wget
wget --spider -r --no-parent \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/suppl/
```

### Using R to see what's available:

```R
library(GEOquery)

# Get list of supplementary files without downloading
supp_info <- getGEOSuppFiles("GSE166259", fetch_files = FALSE)
print(supp_info)
```

## Specific Files for Wilson et al. 2022

For the Wilson paper, you'll want:

1. **Raw sequencing data** (FASTQ files):
   - Available through SRA (linked from GEO)
   - Each GSM sample links to SRA run

2. **Processed data** (likely in supplementary):
   - Molar amount estimates for AML samples
   - Read count matrices
   - Spike-in read counts

## Getting SRA Data (Raw FASTQ files)

```bash
# 1. Install SRA toolkit
# Follow: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

# 2. Get SRA accessions from GEO page
# Look for SRA links in sample records (GSM...)

# 3. Download using prefetch and fastq-dump
prefetch SRR14123456  # Replace with actual SRR number
fastq-dump --split-files --gzip SRR14123456 -O data/raw/HCT116/

# Or use fasterq-dump (faster)
fasterq-dump SRR14123456 -O data/raw/HCT116/
gzip data/raw/HCT116/*.fastq
```

## Troubleshooting

### Issue: "no matches found" error
**Cause**: Shell is trying to expand the `?` as a wildcard
**Solution**: Use quotes around the URL

### Issue: Download times out
**Solution**: Increase timeout
```bash
wget --timeout=600 -O file.tar 'ftp://...'
```

Or in R:
```R
options(timeout = 600)
```

### Issue: Certificate errors with HTTPS
**Solution**: Use FTP instead, or bypass certificate check
```bash
wget --no-check-certificate 'https://...'
```

### Issue: Files are huge
**Solution**: 
1. Download only specific files you need
2. Use processed data instead of raw
3. Download in parts using FTP range requests

## Recommended Workflow

For reproducing the Wilson paper with HCT116 data:

```bash
# 1. Create directory
mkdir -p data/raw/GSE166259

# 2. Download using R (easiest)
R -e 'library(GEOquery); getGEOSuppFiles("GSE166259", baseDir="data/raw/")'

# 3. Check what you got
ls -lh data/raw/GSE166259/

# 4. Extract if needed
cd data/raw/GSE166259
tar -xvf GSE166259_RAW.tar

# 5. Organize by sample type
mkdir -p HCT116 AML
# Move files to appropriate directories based on GSM metadata
```

## Quick Test

Test if FTP access works:

```bash
# Test FTP connection
curl -l ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166259/ | head

# If this works, you can download files
```

---

## Summary: Best Method

For most users, **Method 4 (R/GEOquery)** is recommended because:
- Handles authentication/connections automatically
- Organizes files properly
- Works reliably across platforms
- Part of Bioconductor ecosystem

```R
library(GEOquery)
options(timeout = 600)  # Increase timeout for large files
getGEOSuppFiles("GSE166259", baseDir = "data/raw/")
```
