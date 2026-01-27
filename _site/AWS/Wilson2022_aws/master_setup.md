#!/bin/bash
################################################################################
# Master Setup Script for cfMeDIP-seq Analysis on AWS EC2
# Wilson et al. 2022 Spike-In Normalization Pipeline
#
# Instance: m5.2xlarge (8 vCPUs, 32 GB RAM)
# OS: Ubuntu 22.04 LTS
# 
# Usage: ./master_setup.sh
# 
# This script will:
#   1. Mount data volume
#   2. Install all required software
#   3. Install R packages
#   4. Set up directory structure
#   5. Download reference genomes
#   6. Verify installation
#
# Estimated time: 2-2.5 hours
################################################################################

set -e  # Exit on any error

# Color output for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Log function
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
    exit 1
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

################################################################################
# Configuration Variables
################################################################################

DATA_VOLUME="/dev/nvme1n1"  # Adjust if different (check with lsblk)
DATA_DIR="/data"
HOME_DIR="$HOME"
SCRIPTS_DIR="$HOME/scripts"
PROJECTS_DIR="$HOME/projects"

# Create log directory
mkdir -p $HOME/setup_logs
LOG_FILE="$HOME/setup_logs/master_setup_$(date +'%Y%m%d_%H%M%S').log"

# Redirect all output to log file and console
exec > >(tee -a "$LOG_FILE")
exec 2>&1

################################################################################
# Step 1: Mount Data Volume
################################################################################

mount_data_volume() {
    log "Step 1/7: Mounting data volume"
    
    # Check if data volume exists
    if [ ! -b "$DATA_VOLUME" ]; then
        warn "Data volume $DATA_VOLUME not found. Checking available volumes..."
        lsblk
        warn "Please update DATA_VOLUME variable in this script if needed"
        info "Skipping data volume mount. You can mount manually later."
        return
    fi
    
    # Check if already mounted
    if mountpoint -q "$DATA_DIR"; then
        info "Data volume already mounted at $DATA_DIR"
        return
    fi
    
    # Check if volume is formatted
    if ! sudo file -s "$DATA_VOLUME" | grep -q "ext4"; then
        info "Formatting data volume as ext4..."
        sudo mkfs.ext4 "$DATA_VOLUME"
    fi
    
    # Create mount point
    sudo mkdir -p "$DATA_DIR"
    
    # Mount the volume
    sudo mount "$DATA_VOLUME" "$DATA_DIR"
    sudo chown -R ubuntu:ubuntu "$DATA_DIR"
    
    # Make mount permanent
    if ! grep -q "$DATA_VOLUME" /etc/fstab; then
        echo "$DATA_VOLUME $DATA_DIR ext4 defaults,nofail 0 2" | sudo tee -a /etc/fstab
    fi
    
    # Verify mount
    if mountpoint -q "$DATA_DIR"; then
        log "Data volume successfully mounted at $DATA_DIR"
        df -h "$DATA_DIR"
    else
        error "Failed to mount data volume"
    fi
}

################################################################################
# Step 2: Install System Software
################################################################################

install_system_software() {
    log "Step 2/7: Installing system software and dependencies"
    
    # Update system
    info "Updating system packages..."
    sudo apt update
    sudo apt upgrade -y
    
    # Install essential build tools and libraries
    info "Installing development tools and libraries..."
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
        libjpeg-dev \
        libgit2-dev
    
    # Install R (version 4.3+)
    info "Installing R..."
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
        sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | \
        sudo tee -a /etc/apt/sources.list.d/cran.list
    sudo apt update
    sudo apt install -y r-base r-base-dev
    
    # Install Bowtie2
    info "Installing Bowtie2..."
    sudo apt install -y bowtie2
    
    # Install Samtools
    info "Installing Samtools..."
    sudo apt install -y samtools
    
    # Install Bedtools
    info "Installing Bedtools..."
    sudo apt install -y bedtools
    
    # Install FastQC
    info "Installing FastQC..."
    sudo apt install -y fastqc
    
    # Install SRA Toolkit
    info "Installing SRA Toolkit..."
    cd /tmp
    wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzf sratoolkit.current-ubuntu64.tar.gz
    sudo mv sratoolkit.* /opt/sratoolkit
    
    # Add SRA toolkit to PATH
    if ! grep -q "/opt/sratoolkit/bin" ~/.bashrc; then
        echo 'export PATH=$PATH:/opt/sratoolkit/bin' >> ~/.bashrc
    fi
    export PATH=$PATH:/opt/sratoolkit/bin
    
    # Install additional useful tools
    info "Installing additional tools..."
    sudo apt install -y \
        htop \
        tree \
        ncdu \
        screen \
        tmux
    
    log "System software installation complete"
    
    # Display versions
    info "Installed software versions:"
    echo "  R: $(R --version | head -1)"
    echo "  Bowtie2: $(bowtie2 --version | head -1)"
    echo "  Samtools: $(samtools --version | head -1)"
    echo "  Bedtools: $(bedtools --version | head -1)"
    echo "  FastQC: $(fastqc --version)"
}

################################################################################
# Step 3: Install R Packages
################################################################################

install_r_packages() {
    log "Step 3/7: Installing R packages"
    
    info "This step will take 20-40 minutes..."
    
    cat > /tmp/install_r_packages.R << 'RSCRIPT'
#!/usr/bin/env Rscript
# R package installation for cfMeDIP-seq analysis

cat("=== Installing R Packages ===\n\n")

# Set options
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 4)  # Use multiple cores

# Install BiocManager
cat("[1/5] Installing BiocManager...\n")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)

# Core Bioconductor packages
cat("[2/5] Installing core Bioconductor packages...\n")
BiocManager::install(c(
  "GenomicRanges",
  "GenomicAlignments",
  "Rsamtools",
  "rtracklayer",
  "GenomicFeatures",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Hsapiens.UCSC.hg19"
), ask = FALSE, update = FALSE)

# MEDIPS and methylation tools
cat("[3/5] Installing MEDIPS and methylation tools...\n")
BiocManager::install(c(
  "MEDIPS",
  "edgeR",
  "DNAcopy"
), ask = FALSE, update = FALSE)

# Data manipulation packages
cat("[4/5] Installing tidyverse and data tools...\n")
install.packages(c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "reshape2",
  "cowplot",
  "pheatmap",
  "RColorBrewer"
))

# GEO data tools
cat("[5/5] Installing GEO data tools...\n")
BiocManager::install("GEOquery", ask = FALSE, update = FALSE)

cat("\n=== Package Installation Complete ===\n\n")

# Verify installations
cat("Verifying key packages:\n")
packages <- c("GenomicRanges", "Rsamtools", "MEDIPS", "tidyverse", "GEOquery")
for(pkg in packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  ✓ %s installed\n", pkg))
  } else {
    cat(sprintf("  ✗ %s FAILED\n", pkg))
  }
}
RSCRIPT

    Rscript /tmp/install_r_packages.R
    
    log "R packages installation complete"
}

################################################################################
# Step 4: Set Up Directory Structure
################################################################################

setup_directories() {
    log "Step 4/7: Setting up directory structure"
    
    # Create main directories
    mkdir -p "$DATA_DIR"/{references,raw_data,processed_data,results,scripts}
    mkdir -p "$DATA_DIR"/references/{hg38,hg19,spike_ins,mappability}
    mkdir -p "$DATA_DIR"/raw_data/{fastq,geo_downloads}
    mkdir -p "$DATA_DIR"/processed_data/{bam,bedgraph,counts,qc}
    mkdir -p "$DATA_DIR"/results/{figures,tables}
    
    # Create home directory structure
    mkdir -p "$SCRIPTS_DIR"
    mkdir -p "$PROJECTS_DIR"/wilson_spikein
    
    # Create symbolic link from home
    ln -sf "$DATA_DIR" "$HOME_DIR"/data
    
    # Set permissions
    sudo chown -R ubuntu:ubuntu "$DATA_DIR"
    
    log "Directory structure created"
    tree -L 2 "$DATA_DIR"
}

################################################################################
# Step 5: Download Reference Genome
################################################################################

download_references() {
    log "Step 5/7: Downloading hg38 reference genome"
    
    REF_DIR="$DATA_DIR/references/hg38"
    cd "$REF_DIR"
    
    # Download hg38 reference
    if [ ! -f "hg38.fa.gz" ]; then
        info "Downloading hg38 reference genome (~3 GB)..."
        wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    else
        info "hg38.fa.gz already exists, skipping download"
    fi
    
    # Decompress
    if [ ! -f "hg38.fa" ]; then
        info "Decompressing genome..."
        gunzip -c hg38.fa.gz > hg38.fa
    else
        info "hg38.fa already exists, skipping decompression"
    fi
    
    # Create samtools index
    if [ ! -f "hg38.fa.fai" ]; then
        info "Creating samtools index..."
        samtools faidx hg38.fa
    else
        info "Samtools index already exists"
    fi
    
    # Build Bowtie2 index
    if [ ! -f "hg38.1.bt2" ]; then
        warn "Building Bowtie2 index (~1 hour, 4 GB RAM)..."
        info "This is running in the background. You can monitor with: ps aux | grep bowtie2-build"
        bowtie2-build hg38.fa hg38
    else
        info "Bowtie2 index already exists"
    fi
    
    log "hg38 reference setup complete"
    ls -lh "$REF_DIR"
}

################################################################################
# Step 6: Create Helper Scripts
################################################################################

create_helper_scripts() {
    log "Step 6/7: Creating helper scripts"
    
    # Script to download GEO data
    cat > "$SCRIPTS_DIR/download_geo_data.R" << 'RSCRIPT'
#!/usr/bin/env Rscript
# Download GSE166259 data from Wilson et al. 2022

library(GEOquery)

# Set download directory
setwd("/data/raw_data/geo_downloads")

cat("Downloading GSE166259 metadata...\n")
gse <- getGEO("GSE166259", GSEMatrix = TRUE)

# Get sample information
metadata <- pData(phenoData(gse[[1]]))
write.csv(metadata, "GSE166259_metadata.csv", row.names = FALSE)

cat("Sample information saved to GSE166259_metadata.csv\n")
cat("Number of samples:", nrow(metadata), "\n\n")

# Extract SRA accessions
cat("SRA Run IDs:\n")
if("Run" %in% colnames(metadata)) {
    runs <- metadata$Run
    writeLines(runs, "sra_accessions.txt")
    cat(paste(runs, collapse="\n"))
} else {
    cat("Run column not found. Check metadata columns:\n")
    print(colnames(metadata))
}
RSCRIPT

    # Script to download FASTQ files
    cat > "$SCRIPTS_DIR/download_fastq.sh" << 'BASHSCRIPT'
#!/bin/bash
# Download FASTQ files from SRA

# Usage: ./download_fastq.sh SRR_ACCESSION

if [ $# -eq 0 ]; then
    echo "Usage: $0 SRR_ACCESSION"
    echo "Example: $0 SRR14123456"
    exit 1
fi

SRR=$1
FASTQ_DIR="/data/raw_data/fastq"

cd "$FASTQ_DIR"

echo "Downloading $SRR..."
prefetch "$SRR"

echo "Converting to FASTQ..."
fasterq-dump --split-files "$SRR"

echo "Compressing..."
pigz "$SRR"*.fastq

echo "Complete! Files:"
ls -lh "$SRR"*.fastq.gz
BASHSCRIPT

    # Make scripts executable
    chmod +x "$SCRIPTS_DIR"/*.sh "$SCRIPTS_DIR"/*.R
    
    log "Helper scripts created in $SCRIPTS_DIR"
}

################################################################################
# Step 7: Verify Installation
################################################################################

verify_installation() {
    log "Step 7/7: Verifying installation"
    
    cat > /tmp/verify_setup.R << 'RSCRIPT'
#!/usr/bin/env Rscript
# Verify R packages

packages <- c("GenomicRanges", "Rsamtools", "MEDIPS", 
              "tidyverse", "GEOquery", "edgeR")

cat("Checking R packages:\n")
all_good <- TRUE
for(pkg in packages) {
  status <- require(pkg, character.only=TRUE, quietly=TRUE)
  cat(sprintf("  %s %s\n", pkg, ifelse(status, "✓", "✗")))
  if(!status) all_good <- FALSE
}

if(all_good) {
  cat("\n✓ All R packages installed successfully\n")
  quit(status=0)
} else {
  cat("\n✗ Some R packages failed to install\n")
  quit(status=1)
}
RSCRIPT

    echo ""
    echo "========================================="
    echo "Installation Verification"
    echo "========================================="
    echo ""
    
    # Check software
    echo "Software versions:"
    echo "  R: $(R --version | head -1)"
    echo "  Bowtie2: $(bowtie2 --version | head -1)"
    echo "  Samtools: $(samtools --version | head -1)"
    echo "  Bedtools: $(bedtools --version | head -1)"
    echo "  FastQC: $(fastqc --version)"
    echo ""
    
    # Check R packages
    echo "R packages:"
    Rscript /tmp/verify_setup.R
    echo ""
    
    # Check directories
    echo "Directory structure:"
    tree -L 2 "$DATA_DIR"
    echo ""
    
    # Check disk space
    echo "Disk usage:"
    df -h "$DATA_DIR"
    echo ""
    
    # Check reference files
    echo "Reference files:"
    if [ -f "$DATA_DIR/references/hg38/hg38.fa" ]; then
        echo "  ✓ hg38.fa present ($(stat -c%s "$DATA_DIR/references/hg38/hg38.fa" | numfmt --to=iec))"
    else
        echo "  ✗ hg38.fa missing"
    fi
    
    if [ -f "$DATA_DIR/references/hg38/hg38.1.bt2" ]; then
        echo "  ✓ Bowtie2 index present"
    else
        echo "  ✗ Bowtie2 index missing"
    fi
    echo ""
    
    log "Verification complete!"
}

################################################################################
# Main Execution
################################################################################

main() {
    log "Starting master setup script"
    log "Log file: $LOG_FILE"
    echo ""
    
    # Confirm before proceeding
    echo "This script will set up a complete bioinformatics environment."
    echo "Estimated time: 2-2.5 hours"
    echo ""
    read -p "Continue? (y/n) " -n 1 -r
    echo ""
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Setup cancelled"
        exit 0
    fi
    
    # Record start time
    START_TIME=$(date +%s)
    
    # Run all setup steps
    mount_data_volume
    install_system_software
    install_r_packages
    setup_directories
    download_references
    create_helper_scripts
    verify_installation
    
    # Calculate elapsed time
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    HOURS=$((ELAPSED / 3600))
    MINUTES=$(((ELAPSED % 3600) / 60))
    
    echo ""
    echo "========================================="
    log "Setup complete!"
    echo "========================================="
    echo ""
    echo "Total time: ${HOURS}h ${MINUTES}m"
    echo ""
    echo "Next steps:"
    echo "  1. Review the log file: $LOG_FILE"
    echo "  2. Download GEO data: cd $SCRIPTS_DIR && ./download_geo_data.R"
    echo "  3. Start analysis: cd $PROJECTS_DIR/wilson_spikein"
    echo ""
    echo "Helpful commands:"
    echo "  - Check disk usage: df -h /data"
    echo "  - Monitor processes: htop"
    echo "  - View directory tree: tree -L 3 /data"
    echo ""
}

# Run main function
main "$@"
