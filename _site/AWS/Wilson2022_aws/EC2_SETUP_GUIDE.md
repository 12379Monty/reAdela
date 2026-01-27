# Complete AWS EC2 Setup Guide for cfMeDIP-seq Analysis
## Wilson et al. 2022 Spike-In Normalization Pipeline

**Target System**: AWS EC2 m5.2xlarge (8 vCPUs, 32 GB RAM)  
**Operating System**: Ubuntu 22.04 LTS  
**Estimated Setup Time**: 2-2.5 hours  
**Primary Application**: Cell-free methylated DNA immunoprecipitation sequencing analysis

---

## Table of Contents

1. [Pre-Launch Checklist](#pre-launch-checklist)
2. [Instance Launch Configuration](#instance-launch-configuration)
3. [Initial Connection](#initial-connection)
4. [Automated Setup (Recommended)](#automated-setup-recommended)
5. [Manual Setup (Alternative)](#manual-setup-alternative)
6. [Post-Setup Verification](#post-setup-verification)
7. [Directory Structure](#directory-structure)
8. [Software Inventory](#software-inventory)
9. [Next Steps](#next-steps)

---

## Pre-Launch Checklist

Before launching your EC2 instance, ensure you have:

- [ ] AWS account with appropriate permissions
- [ ] SSH key pair created (or ready to create new one)
- [ ] Security group allowing SSH (port 22) from your IP
- [ ] Optional: Security group allowing HTTP (port 8787) for RStudio Server
- [ ] Budget allocated (~$12-15/day for m5.2xlarge on-demand)

---

## Instance Launch Configuration

### Step 1: Choose AMI
- **AMI**: Ubuntu Server 22.04 LTS (HVM), SSD Volume Type
- **AMI ID**: Search for `ubuntu/images/hvm-ssd/ubuntu-jammy-22.04` in your region

### Step 2: Choose Instance Type
- **Instance Type**: `m5.2xlarge`
- **Specifications**: 8 vCPUs, 32 GB RAM
- **Cost**: ~$0.38/hour (on-demand pricing)

### Step 3: Configure Instance Details
- **Number of instances**: 1
- **Network**: Default VPC
- **Auto-assign Public IP**: Enable
- **IAM role**: (Optional) Create role with S3 read access if storing data in S3

### Step 4: Add Storage

Configure two volumes:

#### Root Volume
- **Size**: 50 GB
- **Volume Type**: gp3
- **Device**: /dev/sda1
- **Delete on Termination**: Yes

#### Data Volume
- **Size**: 300 GB
- **Volume Type**: gp3
- **IOPS**: 3000 (default)
- **Throughput**: 125 MB/s (default)
- **Device**: /dev/sdb (will appear as /dev/nvme1n1)
- **Delete on Termination**: No (preserve data if instance is terminated)

### Step 5: Add Tags
Optional but recommended:
- **Key**: Name | **Value**: cfMeDIP-seq-analysis
- **Key**: Project | **Value**: Wilson-Spikein
- **Key**: Owner | **Value**: [Your Name]

### Step 6: Configure Security Group

Create new or select existing security group:

| Type | Protocol | Port Range | Source | Description |
|------|----------|------------|--------|-------------|
| SSH | TCP | 22 | My IP | SSH access |
| Custom TCP | TCP | 8787 | My IP | RStudio Server (optional) |

**Security Note**: Always restrict SSH to your IP address, not 0.0.0.0/0

### Step 7: Review and Launch
- Review all settings
- Select or create SSH key pair
- **IMPORTANT**: Download and save your private key (.pem file) securely
- Launch instance

---

## Initial Connection

### For macOS/Linux:

```bash
# Set key permissions (first time only)
chmod 400 ~/path/to/your-key.pem

# Connect to instance
ssh -i ~/path/to/your-key.pem ubuntu@YOUR_INSTANCE_PUBLIC_IP
```

### For Windows:

**Option 1: Windows Subsystem for Linux (WSL)**
```bash
chmod 400 /mnt/c/path/to/your-key.pem
ssh -i /mnt/c/path/to/your-key.pem ubuntu@YOUR_INSTANCE_PUBLIC_IP
```

**Option 2: PuTTY**
1. Convert .pem to .ppk using PuTTYgen
2. Use PuTTY to connect with converted key

### Find Your Instance IP

In AWS Console:
1. Navigate to EC2 → Instances
2. Select your instance
3. Copy "Public IPv4 address" from instance details

---

## Automated Setup (Recommended)

### Quick Start

Once connected to your instance, run these commands:

```bash
# Download the master setup script
wget https://raw.githubusercontent.com/[YOUR_REPO]/master_setup.sh

# Make it executable
chmod +x master_setup.sh

# Run the setup
./master_setup.sh
```

The script will:
1. Mount the data volume at `/data`
2. Install all required software (R, Bowtie2, Samtools, etc.)
3. Install R packages (BiocManager, MEDIPS, tidyverse, etc.)
4. Create organized directory structure
5. Download hg38 reference genome
6. Build Bowtie2 indices
7. Create helper scripts
8. Verify installation

**Total Time**: ~2-2.5 hours

### Monitoring Progress

The script outputs to both console and a log file:
```bash
# View the log in real-time (in another terminal)
tail -f ~/setup_logs/master_setup_*.log
```

---

## Manual Setup (Alternative)

If you prefer manual control or the automated script fails, follow these steps:

### 1. Mount Data Volume

```bash
# Check available volumes
lsblk

# Expected output will show /dev/nvme1n1 (or similar)

# Format the volume (WARNING: This erases existing data!)
sudo mkfs.ext4 /dev/nvme1n1

# Create mount point
sudo mkdir /data

# Mount the volume
sudo mount /dev/nvme1n1 /data

# Set ownership
sudo chown ubuntu:ubuntu /data

# Make mount permanent
echo '/dev/nvme1n1 /data ext4 defaults,nofail 0 2' | sudo tee -a /etc/fstab

# Verify
df -h /data
```

### 2. System Update and Essential Tools

```bash
# Update package lists
sudo apt update
sudo apt upgrade -y

# Install build essentials
sudo apt install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    unzip \
    pigz \
    parallel \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev
```

### 3. Install R (Version 4.3+)

```bash
# Add CRAN repository
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# Add R repository
echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | \
    sudo tee -a /etc/apt/sources.list.d/cran.list

# Update and install
sudo apt update
sudo apt install -y r-base r-base-dev

# Verify
R --version
```

### 4. Install Bioinformatics Tools

```bash
# Bowtie2 (alignment)
sudo apt install -y bowtie2
bowtie2 --version

# Samtools (BAM manipulation)
sudo apt install -y samtools
samtools --version

# Bedtools (genomic arithmetic)
sudo apt install -y bedtools
bedtools --version

# FastQC (quality control)
sudo apt install -y fastqc
fastqc --version
```

### 5. Install SRA Toolkit

```bash
# Download and install SRA Toolkit
cd /tmp
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
sudo mv sratoolkit.* /opt/sratoolkit

# Add to PATH
echo 'export PATH=$PATH:/opt/sratoolkit/bin' >> ~/.bashrc
source ~/.bashrc

# Verify
fastq-dump --version
```

### 6. Install R Packages

Create a file `install_packages.R`:

```r
#!/usr/bin/env Rscript

# Set options
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 4)

# Install BiocManager
install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Core Bioconductor packages
BiocManager::install(c(
  "GenomicRanges",
  "GenomicAlignments",
  "Rsamtools",
  "rtracklayer",
  "GenomicFeatures",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Hsapiens.UCSC.hg19"
))

# MEDIPS and methylation tools
BiocManager::install(c("MEDIPS", "edgeR", "DNAcopy"))

# Data manipulation
install.packages(c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "reshape2",
  "cowplot",
  "pheatmap"
))

# GEO data access
BiocManager::install("GEOquery")
```

Run the installation:

```bash
Rscript install_packages.R 2>&1 | tee r_install.log
```

**Time**: 20-40 minutes

### 7. Create Directory Structure

```bash
# Main directories
mkdir -p /data/{references,raw_data,processed_data,results,scripts}

# Reference subdirectories
mkdir -p /data/references/{hg38,hg19,spike_ins,mappability}

# Raw data subdirectories
mkdir -p /data/raw_data/{fastq,geo_downloads}

# Processed data subdirectories
mkdir -p /data/processed_data/{bam,bedgraph,counts,qc}

# Results subdirectories
mkdir -p /data/results/{figures,tables}

# Create symbolic link
ln -s /data ~/data

# Set permissions
sudo chown -R ubuntu:ubuntu /data

# Verify structure
tree -L 2 /data
```

### 8. Download hg38 Reference Genome

```bash
cd /data/references/hg38

# Download reference (3 GB)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Decompress
gunzip -c hg38.fa.gz > hg38.fa

# Create samtools index
samtools faidx hg38.fa

# Build Bowtie2 index (~1 hour)
bowtie2-build hg38.fa hg38
```

---

## Post-Setup Verification

### Verify Software Installations

```bash
# Check all installed software
echo "R: $(R --version | head -1)"
echo "Bowtie2: $(bowtie2 --version | head -1)"
echo "Samtools: $(samtools --version | head -1)"
echo "Bedtools: $(bedtools --version | head -1)"
echo "FastQC: $(fastqc --version)"
```

### Verify R Packages

Create `verify_packages.R`:

```r
#!/usr/bin/env Rscript

packages <- c(
  "GenomicRanges",
  "Rsamtools",
  "MEDIPS",
  "tidyverse",
  "GEOquery",
  "edgeR"
)

cat("Verifying R packages:\n")
for(pkg in packages) {
  status <- require(pkg, character.only = TRUE, quietly = TRUE)
  cat(sprintf("  %s: %s\n", pkg, ifelse(status, "✓", "✗")))
}
```

Run verification:

```bash
Rscript verify_packages.R
```

### Verify Directory Structure

```bash
tree -L 2 /data
```

Expected output:
```
/data
├── references
│   ├── hg38
│   ├── hg19
│   ├── mappability
│   └── spike_ins
├── raw_data
│   ├── fastq
│   └── geo_downloads
├── processed_data
│   ├── bam
│   ├── bedgraph
│   ├── counts
│   └── qc
├── results
│   ├── figures
│   └── tables
└── scripts
```

### Verify Reference Files

```bash
ls -lh /data/references/hg38/

# Expected files:
# - hg38.fa (3.0 GB)
# - hg38.fa.fai (samtools index)
# - hg38.*.bt2 (6 Bowtie2 index files)
```

### Check Disk Usage

```bash
df -h /data

# Should show ~297 GB available (after hg38 download)
```

---

## Directory Structure

### Complete Layout

```
/data/
├── references/               # Reference genomes and annotations
│   ├── hg38/                # Human reference genome (GRCh38)
│   │   ├── hg38.fa          # FASTA sequence
│   │   ├── hg38.fa.fai      # Samtools index
│   │   └── hg38.*.bt2       # Bowtie2 indices
│   ├── hg19/                # Human reference genome (GRCh37)
│   ├── spike_ins/           # Synthetic spike-in controls
│   │   ├── spike_controls.fa
│   │   └── spike_controls.*.bt2
│   └── mappability/         # Mappability tracks
│       └── hg38_mappability.bigWig
│
├── raw_data/                # Raw sequencing data
│   ├── fastq/               # FASTQ files from sequencer
│   │   ├── sample1_R1.fastq.gz
│   │   └── sample1_R2.fastq.gz
│   └── geo_downloads/       # GEO/SRA downloads
│       └── GSE166259/
│
├── processed_data/          # Intermediate processing files
│   ├── bam/                 # Aligned reads
│   │   ├── sample1.bam
│   │   └── sample1.bam.bai
│   ├── bedgraph/            # Coverage tracks
│   ├── counts/              # Read count matrices
│   └── qc/                  # Quality control reports
│       └── fastqc/
│
├── results/                 # Final analysis outputs
│   ├── figures/             # Publication-quality plots
│   │   ├── spike_in_normalization.pdf
│   │   └── methylation_heatmap.pdf
│   └── tables/              # Summary statistics
│       ├── differential_methylation.csv
│       └── spike_in_metrics.csv
│
└── scripts/                 # Analysis scripts
    ├── download_geo_data.R
    ├── download_fastq.sh
    └── align_samples.sh

/home/ubuntu/
├── data -> /data            # Symbolic link for convenience
├── projects/                # Analysis projects
│   └── wilson_spikein/      # Current project
│       ├── README.md
│       ├── analysis.R
│       └── notebooks/
├── scripts/                 # User scripts
└── setup_logs/              # Installation logs
```

---

## Software Inventory

### System Tools

| Software | Version | Purpose |
|----------|---------|---------|
| Ubuntu | 22.04 LTS | Operating system |
| R | 4.3+ | Statistical computing |
| Bowtie2 | 2.4+ | Read alignment |
| Samtools | 1.17+ | BAM file manipulation |
| Bedtools | 2.30+ | Genomic arithmetic |
| FastQC | 0.12+ | Quality control |
| SRA Toolkit | 3.0+ | Download SRA data |
| pigz | 2.6 | Parallel gzip |
| GNU Parallel | 20220522 | Parallel processing |

### R Packages

#### Bioconductor Core
- **GenomicRanges**: Genomic interval operations
- **GenomicAlignments**: Working with aligned reads
- **Rsamtools**: BAM file I/O
- **rtracklayer**: Import/export genomic tracks
- **GenomicFeatures**: Gene annotation
- **BSgenome.Hsapiens.UCSC.hg38**: Human genome sequence
- **BSgenome.Hsapiens.UCSC.hg19**: Human genome sequence (hg19)

#### Methylation Analysis
- **MEDIPS**: Methylated DNA immunoprecipitation sequencing analysis
- **edgeR**: Differential expression/methylation analysis
- **DNAcopy**: Copy number analysis

#### Data Manipulation
- **tidyverse**: Data manipulation suite (dplyr, ggplot2, tidyr, etc.)
- **data.table**: Fast data manipulation
- **reshape2**: Data reshaping

#### Visualization
- **ggplot2**: Publication-quality graphics
- **pheatmap**: Pretty heatmaps
- **cowplot**: Plot composition
- **RColorBrewer**: Color palettes

#### Data Access
- **GEOquery**: Download GEO datasets

---

## Next Steps

### 1. Download Wilson et al. Data

```bash
cd /data/raw_data/geo_downloads
Rscript ~/scripts/download_geo_data.R
```

This will create:
- `GSE166259_metadata.csv`: Sample information
- `sra_accessions.txt`: List of SRA run IDs

### 2. Download Specific Samples

```bash
# Example: Download first sample
SRR=$(head -1 /data/raw_data/geo_downloads/sra_accessions.txt)
~/scripts/download_fastq.sh $SRR
```

### 3. Set Up Spike-In References

```bash
cd /data/references/spike_ins

# Transfer your spike-in FASTA file
# (Use scp from local machine or download from cloud storage)

# Build Bowtie2 index
bowtie2-build spike_controls.fa spike_controls
```

### 4. Begin Analysis

```bash
cd ~/projects/wilson_spikein

# Start R session
R

# Or launch RStudio Server (if installed)
# Navigate to: http://YOUR_INSTANCE_IP:8787
```

---

## Optional: Install RStudio Server

For interactive analysis through a web browser:

```bash
# Install dependencies
sudo apt install -y gdebi-core

# Download RStudio Server
wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.12.1-402-amd64.deb

# Install
sudo gdebi -n rstudio-server-2023.12.1-402-amd64.deb

# Set password for ubuntu user
sudo passwd ubuntu

# Check status
sudo systemctl status rstudio-server
```

Access RStudio at: `http://YOUR_INSTANCE_IP:8787`
- **Username**: ubuntu
- **Password**: [password you set]

**Security Note**: Ensure port 8787 is only accessible from your IP in the security group!

---

## Cost Management

### Pricing Estimate (US East region)

**On-Demand Instance (m5.2xlarge)**:
- Compute: $0.384/hour = $9.22/day
- EBS Storage (350 GB gp3): $0.08/GB/month = $28/month = $0.93/day
- **Total**: ~$10.15/day when running

### Cost Optimization Strategies

1. **Stop (Don't Terminate) When Not In Use**
   ```bash
   # From your local machine
   aws ec2 stop-instances --instance-ids i-1234567890abcdef0
   ```
   - Stopped instances only incur EBS storage costs (~$0.93/day)
   - Can restart with all data and configuration intact

2. **Use Spot Instances**
   - Same performance at 60-90% discount
   - Request spot instance during launch
   - Set max price at on-demand rate ($0.384/hour)
   - **Risk**: Instance may be interrupted if spot price exceeds max

3. **Delete Snapshots**
   - EBS snapshots cost $0.05/GB/month
   - Only keep necessary snapshots

4. **Set Billing Alarms**
   ```bash
   # In AWS Console: CloudWatch → Alarms → Billing
   # Set thresholds: $50, $100, $200
   ```

5. **Data Management**
   - Archive completed analyses to S3 ($0.023/GB/month)
   - Delete intermediate files (BAM files can be regenerated)
   - Compress FASTQ files with pigz

### Typical Usage Patterns

| Scenario | Daily Cost | Monthly Cost |
|----------|-----------|--------------|
| Active analysis (24/7) | $10.15 | $305 |
| Business hours only (8h/day) | $4.00 | $120 |
| Stopped (storage only) | $0.93 | $28 |
| Weekend analysis only | $2.90/week | $12/month |

---

## Troubleshooting

See `EC2_TROUBLESHOOTING_GUIDE.md` for common issues and solutions.

---

## Additional Resources

- **AWS EC2 Documentation**: https://docs.aws.amazon.com/ec2/
- **Bioconductor**: https://www.bioconductor.org/
- **MEDIPS Package**: https://bioconductor.org/packages/MEDIPS/
- **cfMeDIP-seq Protocol**: Shen et al., *Nature Protocols* 2019
- **Wilson et al. 2022**: *Cell Reports Methods* (spike-in normalization)

---

**Document Version**: 1.0  
**Last Updated**: January 2026  
**Author**: Fz Research Documentation  
**License**: MIT
