# EC2 and cfMeDIP-seq Pipeline Troubleshooting Guide

**Target System**: AWS EC2 Ubuntu 22.04  
**Pipeline**: Wilson et al. 2022 Spike-In Normalization

---

## Table of Contents

1. [EC2 Instance Issues](#ec2-instance-issues)
2. [Storage and Mounting Problems](#storage-and-mounting-problems)
3. [Software Installation Issues](#software-installation-issues)
4. [R Package Installation Problems](#r-package-installation-problems)
5. [Data Download Issues](#data-download-issues)
6. [Alignment Problems](#alignment-problems)
7. [Memory and Performance Issues](#memory-and-performance-issues)
8. [Pipeline-Specific Issues](#pipeline-specific-issues)
9. [Network and Connectivity](#network-and-connectivity)
10. [Cost and Billing Issues](#cost-and-billing-issues)

---

## EC2 Instance Issues

### Cannot Connect via SSH

**Symptoms**: 
```
ssh: connect to host ec2-xx-xx-xx-xx.compute.amazonaws.com port 22: Connection timed out
```

**Causes and Solutions**:

1. **Security Group Not Configured**
   ```bash
   # Check security group allows SSH from your IP
   # AWS Console → EC2 → Security Groups → Inbound Rules
   # Should have: Type=SSH, Port=22, Source=Your.IP.Address/32
   ```

2. **Wrong Key File**
   ```bash
   # Verify you're using correct .pem file
   ssh -v -i your-key.pem ubuntu@instance-ip
   # Look for "Offering public key" in verbose output
   ```

3. **Key Permissions Too Open**
   ```bash
   # Fix permissions
   chmod 400 your-key.pem
   ```

4. **Instance Not Running**
   ```bash
   # Check instance state in AWS Console
   # Should show "running" with green circle
   ```

5. **Wrong Username**
   ```bash
   # For Ubuntu AMI, username is 'ubuntu', not 'ec2-user' or 'root'
   ssh -i your-key.pem ubuntu@instance-ip  # Correct
   ssh -i your-key.pem ec2-user@instance-ip  # Wrong for Ubuntu
   ```

### Instance Terminates Immediately After Launch

**Symptoms**: Instance state changes from "pending" to "terminated"

**Causes and Solutions**:

1. **EBS Volume Limit Exceeded**
   - Check AWS Service Quotas for EBS volumes in your region
   - Request limit increase if needed

2. **Insufficient Instance Capacity**
   ```
   Error: InsufficientInstanceCapacity
   ```
   - Try different availability zone
   - Try different instance type
   - Use spot instances as alternative

3. **AMI Not Available in Region**
   - Verify AMI ID is for your current region
   - AMI IDs are region-specific

---

## Storage and Mounting Problems

### Data Volume Not Visible

**Symptoms**: `lsblk` doesn't show expected volume

**Solutions**:

```bash
# List all block devices
lsblk

# Expected output should show:
# nvme0n1 (root volume)
# nvme1n1 (data volume)

# If nvme1n1 is missing:
# 1. Check if volume attached in AWS Console
# 2. Wait 30-60 seconds after instance launch
# 3. Try detaching and reattaching volume

# Reboot instance if volume still not visible
sudo reboot
```

### Mount Failed: "already mounted"

**Symptoms**:
```
mount: /data: /dev/nvme1n1 already mounted
```

**Solution**:

```bash
# Check current mounts
mount | grep nvme1n1

# If already mounted at different location
sudo umount /dev/nvme1n1

# Then mount at /data
sudo mount /dev/nvme1n1 /data
```

### Mount Failed: "wrong fs type"

**Symptoms**:
```
mount: wrong fs type, bad option, bad superblock
```

**Solutions**:

```bash
# Check filesystem type
sudo file -s /dev/nvme1n1

# If output is "data" (no filesystem):
sudo mkfs.ext4 /dev/nvme1n1

# If different filesystem (e.g., xfs):
# Backup data first if volume has existing data!
sudo mkfs.ext4 -F /dev/nvme1n1
```

### Disk Space Full During Setup

**Symptoms**:
```
No space left on device
```

**Solutions**:

```bash
# Check disk usage
df -h

# If /data is full:
# 1. Delete unnecessary files
cd /data
du -sh * | sort -h
# Delete large unnecessary files

# 2. Clean up temp files
rm -rf /tmp/*
rm -rf /data/processed_data/bam/*.sam  # Remove SAM files if BAM exists

# 3. Compress FASTQ files
pigz /data/raw_data/fastq/*.fastq

# If root volume is full:
# Clean apt cache
sudo apt clean
sudo apt autoclean

# Remove old kernels
sudo apt autoremove
```

---

## Software Installation Issues

### apt-get update fails

**Symptoms**:
```
E: Could not get lock /var/lib/apt/lists/lock
```

**Solution**:

```bash
# Another process is using apt
# Wait a few minutes for automatic updates to finish
# Or kill the process:
sudo killall apt apt-get
sudo rm /var/lib/apt/lists/lock
sudo rm /var/cache/apt/archives/lock
sudo rm /var/lib/dpkg/lock*

sudo apt update
```

### Bowtie2 Not Found After Installation

**Symptoms**:
```bash
bowtie2: command not found
```

**Solutions**:

```bash
# Check if installed
dpkg -l | grep bowtie2

# If not installed
sudo apt install -y bowtie2

# If installed but not in PATH
which bowtie2
# Should return: /usr/bin/bowtie2

# Add to PATH if needed
echo 'export PATH=$PATH:/usr/bin' >> ~/.bashrc
source ~/.bashrc
```

### samtools version too old

**Symptoms**:
```
samtools version 1.7 found, but 1.17+ required
```

**Solution**:

```bash
# Remove old version
sudo apt remove samtools

# Install from source
cd /tmp
wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2
tar -xjf samtools-1.19.tar.bz2
cd samtools-1.19
./configure --prefix=/usr/local
make
sudo make install

# Verify
samtools --version
```

---

## R Package Installation Problems

### Cannot install R packages: "permission denied"

**Symptoms**:
```
Warning: cannot remove prior installation of package 'ggplot2'
```

**Solution**:

```r
# Don't use sudo with R package installation
# Run R as regular user, not root

# If libraries are in wrong location:
.libPaths()  # Check library paths

# Set user library
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))
```

### BiocManager installation fails

**Symptoms**:
```
installation of package 'BiocManager' had non-zero exit status
```

**Solutions**:

```bash
# Install system dependencies first
sudo apt install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev

# Then try again in R
install.packages("BiocManager")
```

### GenomicRanges installation fails

**Symptoms**:
```
ERROR: dependencies 'XVector', 'IRanges' are not available
```

**Solution**:

```r
# Install dependencies first
BiocManager::install("BiocGenerics")
BiocManager::install("S4Vectors")
BiocManager::install("IRanges")
BiocManager::install("XVector")

# Then install GenomicRanges
BiocManager::install("GenomicRanges")
```

### BSgenome installation extremely slow

**Symptoms**: Download stalls at large BSgenome packages

**Solution**:

```r
# Use multiple cores
options(Ncpus = 4)

# Try different mirror
options(BioC_mirror = "https://bioconductor.org")

# Or use pre-compiled binary (Ubuntu)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", type = "binary")
```

### MEDIPS installation fails: "limma" required

**Symptoms**:
```
package 'limma' is not available
```

**Solution**:

```r
# Install limma first
BiocManager::install("limma")

# Then MEDIPS
BiocManager::install("MEDIPS")
```

---

## Data Download Issues

### GEOquery fails to download

**Symptoms**:
```
Error in download.file(url, destfile) : cannot open URL
```

**Solutions**:

```r
# Set download timeout
options(timeout = 300)

# Try different method
options(download.file.method = "wget")

# Or use curl
options(download.file.method = "curl")

# Retry with GEOquery
library(GEOquery)
gse <- getGEO("GSE166259")
```

### SRA Toolkit: "prefetch" fails

**Symptoms**:
```
failed to resolve accession 'SRR14123456'
```

**Solutions**:

```bash
# Configure SRA Toolkit
vdb-config --interactive
# Enable "Remote Access" and set cache location

# Or use fastq-dump directly
fastq-dump --split-files --gzip SRR14123456

# Alternative: Download from ENA (faster)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/056/SRR14123456/SRR14123456_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR141/056/SRR14123456/SRR14123456_2.fastq.gz
```

### wget fails: "Connection refused"

**Symptoms**:
```
Connecting to hgdownload.soe.ucsc.edu... failed: Connection refused
```

**Solutions**:

```bash
# Check network connectivity
ping 8.8.8.8

# If ping works, try curl instead
curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Use --retry option
wget --retry-connrefused --waitretry=1 --timeout=20 --tries=10 \
    http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

---

## Alignment Problems

### Bowtie2 index not found

**Symptoms**:
```
Could not locate a Bowtie index corresponding to basename "hg38"
```

**Solutions**:

```bash
# Check index files exist
ls -lh /data/references/hg38/hg38*.bt2

# Should see 6 files:
# hg38.1.bt2, hg38.2.bt2, hg38.3.bt2, hg38.4.bt2, hg38.rev.1.bt2, hg38.rev.2.bt2

# If missing, rebuild index
cd /data/references/hg38
bowtie2-build hg38.fa hg38

# In alignment command, use full path to index
bowtie2 -x /data/references/hg38/hg38 ...
```

### Bowtie2 runs out of memory

**Symptoms**:
```
Error: Ran out of memory
terminate called after throwing an instance of 'std::bad_alloc'
```

**Solutions**:

```bash
# Monitor memory usage
htop

# Reduce number of threads
bowtie2 -x hg38 -1 R1.fq -2 R2.fq --threads 4  # Instead of 8

# Upgrade to memory-optimized instance
# r5.2xlarge (64 GB RAM) instead of m5.2xlarge (32 GB)

# Process samples sequentially instead of parallel
for sample in sample1 sample2 sample3; do
    bowtie2 -x hg38 -1 ${sample}_R1.fq -2 ${sample}_R2.fq --threads 8
done
```

### Low alignment rate

**Symptoms**: 
```
35.67% overall alignment rate
```
Expected: >80%

**Causes and Solutions**:

1. **Wrong reference genome**
   ```bash
   # Verify genome version matches data
   # If data is hg19, use hg19 reference, not hg38
   ```

2. **Adapter contamination**
   ```bash
   # Trim adapters first
   # Install cutadapt
   pip install cutadapt --break-system-packages
   
   # Trim adapters
   cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
       -o trimmed_R1.fq -p trimmed_R2.fq \
       input_R1.fq input_R2.fq
   ```

3. **Poor quality reads**
   ```bash
   # Check FastQC report
   # Filter low-quality reads with fastp
   ```

### SAM/BAM files corrupted

**Symptoms**:
```
[E::sam_parse1] missing SAM header
[main_samview] truncated file
```

**Solutions**:

```bash
# Check file integrity
samtools quickcheck file.bam
# If error, file is corrupted

# Re-align sample
# Make sure alignment completes successfully

# Check disk space during alignment
df -h  # If full, alignment may have failed silently
```

---

## Memory and Performance Issues

### R session crashes: "Cannot allocate vector"

**Symptoms**:
```r
Error: cannot allocate vector of size 8.5 Gb
```

**Solutions**:

```r
# Check available memory
gc()

# Process chromosomes separately
for(chr in paste0("chr", 1:22)) {
    # Process one chromosome at a time
    chr_data <- count_matrix[seqnames == chr, ]
    # Analyze
    # Save results
    # Clear memory
    rm(chr_data)
    gc()
}

# Use data.table for large datasets
library(data.table)
dt <- fread("large_file.csv")  # Faster than read.csv

# Increase instance RAM
# Upgrade to r5.2xlarge (64 GB) or r5.4xlarge (128 GB)
```

### Script killed without error message

**Symptoms**: Process terminates with "Killed"

**Cause**: Out of Memory (OOM) killer

**Solutions**:

```bash
# Check system logs
dmesg | grep -i "killed process"

# Monitor memory during execution
# In another terminal:
watch -n 1 free -h

# Solutions:
# 1. Reduce memory usage in script
# 2. Process fewer samples at once
# 3. Upgrade instance type
# 4. Add swap space (temporary fix)

# Add swap (4 GB)
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
# Make permanent
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

### Bowtie2 building index very slow

**Symptoms**: `bowtie2-build` running for 2+ hours

**Expected**: ~60 minutes on m5.2xlarge

**Solutions**:

```bash
# Check if it's actually running
ps aux | grep bowtie2-build

# Monitor progress
# Bowtie2-build doesn't show progress, but you can check memory usage
htop

# If stuck:
# Kill and restart
killall bowtie2-build

# Run with screen so you can detach
screen -S bowtie_build
bowtie2-build hg38.fa hg38
# Detach: Ctrl+A, then D
# Reattach: screen -r bowtie_build

# Or use nohup
nohup bowtie2-build hg38.fa hg38 > bowtie_build.log 2>&1 &
```

---

## Pipeline-Specific Issues

### No spike-in reads detected

**Symptoms**: All samples show 0 spike-in aligned reads

**Causes and Solutions**:

1. **Spike-ins not added to samples**
   - Verify spike-ins were added during library prep
   - Check lab notebook/protocol

2. **Wrong spike-in reference**
   ```bash
   # Verify spike-in sequences
   head /data/references/spike_ins/spike_controls.fa
   
   # Should match sequences from Wilson et al. 2022 Supplementary Data
   ```

3. **Alignment parameters too stringent**
   ```bash
   # Try more lenient alignment
   bowtie2 --very-sensitive  # Instead of --very-sensitive-local
   ```

### Spike-in alignment rate too high (>10%)

**Symptoms**: Spike-in reads are >10% of total

**Likely Cause**: Contamination or incorrect reference

**Solutions**:

```bash
# Check spike-in reference for contamination
# Should only contain synthetic sequences, not genomic DNA

# Verify FASTQ files
# Check for unexpected adapter sequences or contamination

# Re-download spike-in sequences from original source
```

### Normalization factors too variable

**Symptoms**: CV of normalization factors >50%

**Indicates**: Poor batch consistency

**Solutions**:

1. **Review sample processing**
   - Check IP protocol consistency
   - Verify antibody lot numbers
   - Review library prep dates

2. **Check for outliers**
   ```r
   # Identify outlier samples
   outliers <- norm_factors %>%
       filter(norm_factor > median(norm_factor) + 2*sd(norm_factor) |
              norm_factor < median(norm_factor) - 2*sd(norm_factor))
   
   # Remove outliers and re-analyze
   ```

3. **Apply additional normalization**
   ```r
   # Use quantile normalization in addition to spike-in
   library(preprocessCore)
   normalized_matrix <- normalize.quantiles(count_matrix)
   ```

### Coverage tracks show no signal

**Symptoms**: bedGraph files contain all zeros

**Causes and Solutions**:

1. **Wrong normalization factor**
   ```bash
   # Check normalization factors
   cat /data/results/tables/spikein_normalization_factors.csv
   
   # Scale factors should be 0.5-2.0
   # If scale_factor = 0, division by zero occurred
   ```

2. **BAM file issues**
   ```bash
   # Check BAM has reads
   samtools view -c sample.bam
   # Should return >1000000
   
   # Check BAM is sorted
   samtools view -H sample.bam | grep "SO:"
   # Should show SO:coordinate
   ```

3. **bedtools genomecov parameters**
   ```bash
   # Use correct parameters
   bedtools genomecov -ibam sample.bam -bg -scale 1.5
   # -bg: bedGraph format
   # -scale: normalization factor
   ```

---

## Network and Connectivity

### Cannot access external resources

**Symptoms**:
```
curl: (6) Could not resolve host
```

**Solutions**:

```bash
# Check DNS resolution
nslookup google.com

# If DNS fails, add public DNS servers
echo "nameserver 8.8.8.8" | sudo tee -a /etc/resolv.conf
echo "nameserver 8.8.4.4" | sudo tee -a /etc/resolv.conf

# Check security group allows outbound traffic
# AWS Console → EC2 → Security Groups → Outbound Rules
# Should allow all outbound traffic (default)
```

### RStudio Server not accessible

**Symptoms**: Cannot connect to http://instance-ip:8787

**Solutions**:

1. **Security group not configured**
   ```
   AWS Console → EC2 → Security Groups → Inbound Rules
   Add: Type=Custom TCP, Port=8787, Source=Your.IP/32
   ```

2. **RStudio Server not running**
   ```bash
   sudo systemctl status rstudio-server
   
   # If not running:
   sudo systemctl start rstudio-server
   sudo systemctl enable rstudio-server
   ```

3. **Firewall blocking**
   ```bash
   # Check UFW status
   sudo ufw status
   
   # If active, allow port 8787
   sudo ufw allow 8787/tcp
   ```

---

## Cost and Billing Issues

### Unexpected charges

**Common Causes**:

1. **Forgot to stop instance**
   - Running 24/7: $9.22/day
   - **Solution**: Stop when not in use

2. **EBS volumes not deleted**
   - Detached volumes still incur charges
   - **Check**: EC2 → Volumes → Filter: "available"
   - **Delete unused volumes**

3. **Snapshots accumulating**
   - $0.05/GB/month
   - **Check**: EC2 → Snapshots
   - **Delete old snapshots**

4. **Data transfer costs**
   - Downloads from EC2 to internet
   - **Solution**: Process data on EC2, only download results

**Prevention**:

```bash
# Set up billing alerts
# AWS Console → Billing → Budgets
# Create budget: $100/month with 80% alert

# Tag resources for tracking
aws ec2 create-tags --resources i-xxxxx \
    --tags Key=Project,Value=cfMeDIP-seq Key=Owner,Value=Fz
```

### How to stop instance properly

```bash
# From local machine (using AWS CLI)
aws ec2 stop-instances --instance-ids i-xxxxx

# Or from AWS Console
# EC2 → Instances → Select instance → Instance State → Stop

# To restart later
aws ec2 start-instances --instance-ids i-xxxxx
```

---

## Getting Help

### Collecting Information for Support

When reporting issues, collect:

```bash
# System information
uname -a
cat /etc/os-release

# Installed software versions
R --version
bowtie2 --version
samtools --version

# Memory and disk
free -h
df -h

# Error messages
# Copy complete error text

# Relevant log files
cat ~/setup_logs/master_setup_*.log
```

### Useful Commands for Debugging

```bash
# Check process status
ps aux | grep bowtie2
top -u ubuntu

# Check system logs
sudo journalctl -xe
dmesg | tail -50

# Check R session info
R -e "sessionInfo()"

# Network diagnostics
ping 8.8.8.8
curl -I https://cloud.r-project.org

# Disk I/O stats
iostat -x 1 10
```

---

## Quick Reference: Common Fixes

| Problem | Quick Fix |
|---------|-----------|
| Can't SSH | `chmod 400 key.pem` |
| Volume not mounted | `sudo mount /dev/nvme1n1 /data` |
| Command not found | `source ~/.bashrc` |
| Out of disk space | `rm /data/processed_data/bam/*.sam` |
| R package won't install | `sudo apt install -y lib*-dev` |
| Process killed | Upgrade to larger instance |
| Slow download | Use `curl` or try different mirror |
| Bowtie2 index missing | Check path: `/data/references/hg38/hg38` |
| High AWS bill | Stop instance when not in use |

---

## Emergency Procedures

### Instance Unresponsive

```bash
# 1. Try to connect via AWS Console
# EC2 → Instances → Connect → EC2 Instance Connect

# 2. Reboot instance
aws ec2 reboot-instances --instance-ids i-xxxxx

# 3. If still unresponsive, force stop and restart
aws ec2 stop-instances --instance-ids i-xxxxx --force
aws ec2 start-instances --instance-ids i-xxxxx
```

### Data Recovery

```bash
# Create EBS snapshot before making major changes
aws ec2 create-snapshot --volume-id vol-xxxxx \
    --description "Before pipeline run"

# If needed, restore from snapshot:
# 1. Create volume from snapshot
# 2. Attach to instance
# 3. Mount and recover data
```

---

**Document Version**: 1.0  
**Last Updated**: January 2026  
**Maintainer**: Fz Research Documentation
