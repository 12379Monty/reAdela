# EC2 cfMeDIP-seq Quick Reference Guide

**Instance**: m5.2xlarge | **OS**: Ubuntu 22.04 | **Pipeline**: Wilson et al. 2022

---

## Essential Commands

### Instance Management

```bash
# Connect to instance
ssh -i your-key.pem ubuntu@instance-ip

# Check instance status (from local machine)
aws ec2 describe-instances --instance-ids i-xxxxx

# Stop instance (save costs when not in use)
aws ec2 stop-instances --instance-ids i-xxxxx

# Start stopped instance
aws ec2 start-instances --instance-ids i-xxxxx

# Reboot instance
aws ec2 reboot-instances --instance-ids i-xxxxx
```

### Session Management

```bash
# Use screen for long-running processes
screen -S analysis          # Start new session
# Run your commands
Ctrl+A then D              # Detach
screen -r analysis         # Reattach
screen -ls                 # List sessions

# Alternative: tmux
tmux new -s analysis
tmux attach -t analysis
tmux ls
```

---

## Directory Structure

```
/data/                              # Main data directory (300 GB)
├── references/
│   ├── hg38/                       # Human reference genome
│   │   ├── hg38.fa                 # FASTA sequence (3 GB)
│   │   ├── hg38.fa.fai             # Samtools index
│   │   └── hg38.*.bt2              # Bowtie2 indices (6 files)
│   ├── spike_ins/                  # Synthetic controls
│   │   ├── spike_controls.fa
│   │   └── spike_controls.*.bt2
│   └── mappability/
│       └── hg38_mappability.bigWig
│
├── raw_data/
│   ├── fastq/                      # Raw sequencing files
│   │   ├── sample1_R1.fastq.gz
│   │   └── sample1_R2.fastq.gz
│   └── geo_downloads/              # GEO/SRA data
│
├── processed_data/
│   ├── spikein_alignment/          # Spike-in aligned BAMs
│   ├── bam/                        # Genome-aligned BAMs
│   ├── bedgraph/                   # Coverage tracks
│   ├── counts/                     # Count matrices
│   └── qc/                         # Quality reports
│
├── results/
│   ├── figures/                    # Publication plots
│   └── tables/                     # Results tables
│
└── scripts/                        # Analysis scripts

~/                                  # Home directory
├── data -> /data                   # Symbolic link
├── projects/
│   └── wilson_spikein/             # Current project
└── scripts/
    ├── download_geo_data.R
    └── download_fastq.sh
```

---

## Pipeline Workflow

### Quick Start

```bash
# 1. Quality control
./01_fastqc.sh

# 2. Align to spike-ins
./02_align_spikein.sh

# 3. Align to genome
./03_align_genome.sh

# 4. Calculate normalization
Rscript 04_calculate_normalization.R

# 5. Generate coverage
./05_generate_coverage.sh

# 6. Quantify windows
Rscript 06_window_quantification.R

# 7. Differential analysis
Rscript 07_medips_analysis.R
```

### Process Single Sample

```bash
SAMPLE="sample1"
FASTQ_DIR="/data/raw_data/fastq"
SPIKEIN_REF="/data/references/spike_ins/spike_controls"
GENOME_REF="/data/references/hg38/hg38"

# Align to spike-ins
bowtie2 -x "$SPIKEIN_REF" \
    -1 "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
    -2 "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" \
    --threads 8 --very-sensitive \
    --un-conc-gz "/data/processed_data/spikein_alignment/${SAMPLE}_unmapped.fastq.gz" \
    -S "/data/processed_data/spikein_alignment/${SAMPLE}_spikein.sam"

# Align to genome
bowtie2 -x "$GENOME_REF" \
    -1 "/data/processed_data/spikein_alignment/${SAMPLE}_unmapped.fastq.1.gz" \
    -2 "/data/processed_data/spikein_alignment/${SAMPLE}_unmapped.fastq.2.gz" \
    --threads 8 --very-sensitive \
    -S "/data/processed_data/bam/${SAMPLE}_genome.sam"

# Convert to BAM
samtools view -@ 8 -bS -q 10 "/data/processed_data/bam/${SAMPLE}_genome.sam" | \
    samtools sort -@ 8 -o "/data/processed_data/bam/${SAMPLE}.bam"

samtools index "/data/processed_data/bam/${SAMPLE}.bam"
```

---

## Data Management

### Download GEO Data

```r
# In R
library(GEOquery)
setwd("/data/raw_data/geo_downloads")

# Download metadata
gse <- getGEO("GSE166259")
metadata <- pData(phenoData(gse[[1]]))
write.csv(metadata, "GSE166259_metadata.csv")
```

```bash
# Download FASTQ from SRA
SRR="SRR14123456"
cd /data/raw_data/fastq

# Method 1: fastq-dump
fastq-dump --split-files --gzip "$SRR"

# Method 2: fasterq-dump + pigz (faster)
fasterq-dump --split-files "$SRR"
pigz ${SRR}*.fastq

# Method 3: prefetch + fasterq-dump
prefetch "$SRR"
fasterq-dump --split-files "$SRR"
pigz ${SRR}*.fastq
```

### Check File Integrity

```bash
# Verify FASTQ files
gunzip -t file.fastq.gz

# Count reads
zcat file.fastq.gz | echo $((`wc -l`/4))

# Verify BAM files
samtools quickcheck file.bam
samtools view -c file.bam  # Count reads

# Check bedGraph files
head -20 file.bedGraph
```

### Compress/Decompress

```bash
# Compress with pigz (parallel gzip)
pigz file.fastq         # Creates file.fastq.gz
pigz -d file.fastq.gz   # Decompress

# Compress BAM (already compressed)
# No additional compression needed

# Compress results
tar -czf results_20260117.tar.gz /data/results/
```

---

## Monitoring and Diagnostics

### System Resources

```bash
# CPU and memory
htop                    # Interactive
top                     # Basic

# Disk usage
df -h                   # All filesystems
df -h /data             # Just data volume
du -sh /data/*          # Size of each directory
ncdu /data              # Interactive disk usage

# Disk I/O
iostat -x 1 10          # 10 samples, 1 second apart

# Network
iftop                   # Network bandwidth (needs install)
nethogs                 # Per-process bandwidth
```

### Process Monitoring

```bash
# List processes
ps aux | grep bowtie2
ps aux | grep R

# Process tree
pstree -p

# Kill processes
killall bowtie2         # Kill all bowtie2 processes
kill -9 PID             # Force kill specific process

# Background jobs
jobs                    # List background jobs
fg %1                   # Bring job 1 to foreground
bg %1                   # Send job 1 to background
```

### Log Files

```bash
# Setup logs
ls -lht ~/setup_logs/

# View latest setup log
tail -f ~/setup_logs/master_setup_*.log

# Pipeline logs
tail -f /data/processed_data/bam/*_alignment.log

# System logs
sudo journalctl -xe     # Recent errors
dmesg | tail -50        # Kernel messages
```

---

## R Quick Reference

### Launch R

```bash
# Command line R
R

# RStudio Server (if installed)
# Open browser: http://instance-ip:8787
# Login: ubuntu / [your password]
```

### Common R Commands

```r
# Check installed packages
installed.packages()[, "Package"]

# Check library paths
.libPaths()

# Load libraries
library(GenomicRanges)
library(tidyverse)
library(MEDIPS)

# Session information
sessionInfo()

# Memory management
gc()                    # Garbage collection
object.size(x)         # Object size
rm(list = ls())        # Clear workspace

# Set working directory
setwd("/data/projects/wilson_spikein")
getwd()

# Quit R
q()                    # Prompts to save workspace
q(save = "no")         # Don't save
```

### Debug R Scripts

```r
# Run script with error handling
tryCatch({
    source("analysis.R")
}, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    traceback()
})

# Run line by line in R
source("script.R", echo = TRUE)

# Check for errors
warnings()
```

---

## File Transfer

### Upload to EC2

```bash
# From local machine

# Single file
scp -i your-key.pem local_file.txt ubuntu@instance-ip:/data/raw_data/

# Directory
scp -i your-key.pem -r local_directory/ ubuntu@instance-ip:/data/

# With compression
tar czf - local_directory/ | \
    ssh -i your-key.pem ubuntu@instance-ip "cd /data && tar xzf -"
```

### Download from EC2

```bash
# From local machine

# Single file
scp -i your-key.pem ubuntu@instance-ip:/data/results/tables/results.csv .

# Directory
scp -i your-key.pem -r ubuntu@instance-ip:/data/results/ .

# With compression
ssh -i your-key.pem ubuntu@instance-ip "cd /data && tar czf - results/" | \
    tar xzf - -C local_directory/
```

### Using rsync (Resumable)

```bash
# Upload (from local)
rsync -avz -e "ssh -i your-key.pem" \
    local_data/ ubuntu@instance-ip:/data/raw_data/

# Download (from local)
rsync -avz -e "ssh -i your-key.pem" \
    ubuntu@instance-ip:/data/results/ local_results/

# Dry run (test without copying)
rsync -avzn -e "ssh -i your-key.pem" \
    local_data/ ubuntu@instance-ip:/data/
```

---

## Troubleshooting Quick Fixes

| Problem | Quick Fix |
|---------|-----------|
| **Can't SSH** | `chmod 400 key.pem` |
| **Command not found** | `source ~/.bashrc` or `export PATH=$PATH:/usr/bin` |
| **Volume not mounted** | `sudo mount /dev/nvme1n1 /data` |
| **Out of space** | `rm /data/processed_data/bam/*.sam` |
| **Out of memory** | Reduce threads: `--threads 4` instead of 8 |
| **Bowtie2 index error** | Use full path: `-x /data/references/hg38/hg38` |
| **R package error** | Install dependencies: `sudo apt install -y libcurl4-openssl-dev` |
| **Process killed** | Check memory: `dmesg \| grep -i killed` |
| **Slow alignment** | Use screen/tmux and run overnight |

---

## Performance Tips

### Parallel Processing

```bash
# GNU Parallel for batch processing
# Create sample list
ls /data/raw_data/fastq/*_R1.fastq.gz | \
    xargs -n 1 basename | sed 's/_R1.fastq.gz//' > samples.txt

# Run in parallel (4 at a time)
parallel -j 4 './process_sample.sh {}' :::: samples.txt

# With log file
parallel -j 4 --joblog parallel.log \
    './process_sample.sh {}' :::: samples.txt
```

### Memory Optimization

```bash
# Monitor memory before running
free -h

# For large files, process in chunks
# In R:
library(data.table)
fread("large_file.csv", nrows = 1000000)  # First 1M rows

# Use pigz instead of gzip (parallel)
pigz -p 8 large_file.fastq  # Use 8 cores
```

### Speed Up Alignment

```bash
# Use more threads
bowtie2 --threads 8 ...

# Use faster alignment mode (less sensitive)
bowtie2 --fast ...  # Instead of --very-sensitive

# Skip duplicate marking if not needed
# (Duplicates are rare in cfDNA)
```

---

## Backup and Snapshots

### Create EBS Snapshot

```bash
# From local machine with AWS CLI
aws ec2 create-snapshot \
    --volume-id vol-xxxxx \
    --description "Before analysis - $(date +%Y%m%d)"

# List snapshots
aws ec2 describe-snapshots --owner-ids self
```

### Save Results to S3

```bash
# Install AWS CLI (if not already)
sudo apt install -y awscli

# Configure (one time)
aws configure

# Upload results
aws s3 cp /data/results/ s3://my-bucket/wilson-analysis/ --recursive

# Download later
aws s3 cp s3://my-bucket/wilson-analysis/ /data/results/ --recursive
```

---

## Cost Optimization

### Current Costs (US East)

- **Running m5.2xlarge**: $0.384/hour = $9.22/day
- **Storage (350 GB gp3)**: $28/month = $0.93/day
- **Stopped instance**: $0.93/day (storage only)

### Save Money

```bash
# 1. Stop when not in use
aws ec2 stop-instances --instance-ids i-xxxxx

# 2. Delete unnecessary files
rm /data/processed_data/bam/*.sam  # Keep only BAM
rm /data/raw_data/fastq/*.fastq    # Keep only .gz

# 3. Compress results before downloading
tar czf results.tar.gz /data/results/

# 4. Use spot instances (60-90% discount)
# Select "Spot instance" during launch

# 5. Delete old snapshots
aws ec2 delete-snapshot --snapshot-id snap-xxxxx
```

---

## Essential Software Versions

| Software | Version | Check Command |
|----------|---------|---------------|
| Ubuntu | 22.04 | `lsb_release -a` |
| R | ≥4.3 | `R --version` |
| Bowtie2 | ≥2.4 | `bowtie2 --version` |
| Samtools | ≥1.17 | `samtools --version` |
| Bedtools | ≥2.30 | `bedtools --version` |
| FastQC | ≥0.12 | `fastqc --version` |

---

## Important File Paths

```bash
# Reference files
/data/references/hg38/hg38.fa
/data/references/hg38/hg38.*.bt2
/data/references/spike_ins/spike_controls.fa
/data/references/mappability/hg38_mappability.bigWig

# Configuration files
~/.bashrc                           # Bash configuration
~/.Rprofile                         # R configuration
/etc/fstab                          # Mount configuration

# Log files
~/setup_logs/                       # Setup logs
/data/processed_data/bam/*_alignment.log  # Alignment logs

# Results
/data/results/tables/spikein_normalization_factors.csv
/data/results/tables/differential_methylation_results.csv
/data/results/figures/
```

---

## Emergency Contacts and Resources

### Documentation

- **Complete Setup Guide**: `~/EC2_SETUP_GUIDE.md`
- **Pipeline Details**: `~/WILSON_PIPELINE.md`
- **Troubleshooting**: `~/EC2_TROUBLESHOOTING_GUIDE.md`
- **This Quick Reference**: `~/QUICK_REFERENCE.md`

### External Resources

- **AWS EC2 Docs**: https://docs.aws.amazon.com/ec2/
- **Bioconductor**: https://bioconductor.org/
- **MEDIPS Manual**: https://bioconductor.org/packages/MEDIPS/
- **Wilson et al. 2022**: https://doi.org/10.1016/j.crmeth.2022.100366
- **cfMeDIP-seq Protocol**: https://doi.org/10.1038/s41596-019-0202-2

### Support Commands

```bash
# Generate system report for troubleshooting
cat > system_report.txt << EOF
=== System Information ===
$(uname -a)
$(lsb_release -a)

=== Software Versions ===
R: $(R --version | head -1)
Bowtie2: $(bowtie2 --version | head -1)
Samtools: $(samtools --version | head -1)

=== Disk Usage ===
$(df -h)

=== Memory ===
$(free -h)

=== Running Processes ===
$(ps aux | grep -E 'bowtie2|R|samtools')
EOF

cat system_report.txt
```

---

## One-Liners

```bash
# Count total reads in FASTQ
zcat file.fastq.gz | echo $((`wc -l`/4))

# Count aligned reads in BAM
samtools view -c -F 4 file.bam

# Top 10 largest files
du -ah /data | sort -rh | head -10

# Find all BAM files
find /data -name "*.bam" -type f

# Check all running processes
ps aux | grep -v grep | grep -E 'bowtie2|samtools|R'

# Memory usage by process
ps aux --sort=-%mem | head -10

# Count files in directory
find /data/raw_data/fastq -type f | wc -l

# Replace text in all files
find . -name "*.R" -exec sed -i 's/old_text/new_text/g' {} +
```

---

**Quick Reference Version**: 1.0  
**Last Updated**: January 2026  
**For**: Wilson et al. 2022 Spike-In Normalization Pipeline
