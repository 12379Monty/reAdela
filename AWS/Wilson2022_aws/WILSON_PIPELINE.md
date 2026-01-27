# Wilson et al. 2022 Spike-In Normalization Pipeline
## cfMeDIP-seq Analysis Workflow

**Reference**: Wilson et al., *Cell Reports Methods* 2022  
**DOI**: 10.1016/j.crmeth.2022.100366  
**GEO Accession**: GSE166259

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Prerequisites](#prerequisites)
3. [Workflow Steps](#workflow-steps)
4. [Scripts](#scripts)
5. [Expected Outputs](#expected-outputs)
6. [Quality Control Metrics](#quality-control-metrics)

---

## Pipeline Overview

### Conceptual Workflow

```
FASTQ files (cell-free DNA)
    ↓
Quality Control (FastQC)
    ↓
Align to spike-in reference (Bowtie2)
    ↓
Extract unmapped reads
    ↓
Align to human genome (Bowtie2)
    ↓
Convert to BAM, sort, index
    ↓
Calculate spike-in normalization factors
    ↓
Generate normalized coverage tracks
    ↓
Binning and feature counting
    ↓
Statistical analysis & visualization
```

### Key Innovation

The Wilson et al. method uses **synthetic spike-in controls** with known methylation patterns to:
- Normalize for immunoprecipitation efficiency variation
- Correct for batch effects
- Enable quantitative comparison across samples
- Improve reproducibility

---

## Prerequisites

### Required Files

1. **Reference Genomes**
   - `/data/references/hg38/hg38.fa` (+ Bowtie2 indices)
   - `/data/references/spike_ins/spike_controls.fa` (+ Bowtie2 indices)

2. **Mappability Track**
   - `/data/references/mappability/hg38_mappability.bigWig`

3. **Sample Data**
   - Raw FASTQ files in `/data/raw_data/fastq/`
   - Sample metadata sheet

### Software Requirements

All installed by master setup script:
- Bowtie2 (≥2.4.0)
- Samtools (≥1.17)
- Bedtools (≥2.30)
- R (≥4.3)
- R packages: GenomicRanges, Rsamtools, MEDIPS, tidyverse

---

## Workflow Steps

### Step 1: Quality Control

**Script**: `01_fastqc.sh`

```bash
#!/bin/bash
# Quality control of raw FASTQ files using FastQC

FASTQ_DIR="/data/raw_data/fastq"
QC_DIR="/data/processed_data/qc/fastqc"

mkdir -p "$QC_DIR"

# Run FastQC on all FASTQ files
for fastq in "$FASTQ_DIR"/*.fastq.gz; do
    echo "Processing: $(basename $fastq)"
    fastqc -o "$QC_DIR" -t 4 "$fastq"
done

# Generate MultiQC report (if installed)
if command -v multiqc &> /dev/null; then
    multiqc "$QC_DIR" -o "$QC_DIR" -n fastqc_report
fi

echo "FastQC complete. Reports in: $QC_DIR"
```

**Expected Output**:
- Individual HTML reports for each FASTQ file
- Summary statistics on read quality, adapter content, duplication

---

### Step 2: Align to Spike-In Controls

**Script**: `02_align_spikein.sh`

```bash
#!/bin/bash
# Align reads to synthetic spike-in controls
# Unmapped reads will be used for downstream genome alignment

set -e

FASTQ_DIR="/data/raw_data/fastq"
SPIKEIN_REF="/data/references/spike_ins/spike_controls"
OUTPUT_DIR="/data/processed_data/spikein_alignment"
THREADS=8

mkdir -p "$OUTPUT_DIR"

# Function to align paired-end reads
align_spikein() {
    local sample=$1
    local r1="${FASTQ_DIR}/${sample}_R1.fastq.gz"
    local r2="${FASTQ_DIR}/${sample}_R2.fastq.gz"
    
    echo "=== Processing sample: $sample ==="
    
    # Align to spike-in controls, save unmapped reads
    bowtie2 \
        -x "$SPIKEIN_REF" \
        -1 "$r1" \
        -2 "$r2" \
        --threads "$THREADS" \
        --very-sensitive \
        --no-unal \
        --un-conc-gz "${OUTPUT_DIR}/${sample}_unmapped.fastq.gz" \
        -S "${OUTPUT_DIR}/${sample}_spikein.sam" \
        2> "${OUTPUT_DIR}/${sample}_spikein_alignment.log"
    
    # Convert to BAM, sort, and index
    samtools view -@ "$THREADS" -bS "${OUTPUT_DIR}/${sample}_spikein.sam" | \
        samtools sort -@ "$THREADS" -o "${OUTPUT_DIR}/${sample}_spikein.bam"
    
    samtools index "${OUTPUT_DIR}/${sample}_spikein.bam"
    
    # Remove SAM file to save space
    rm "${OUTPUT_DIR}/${sample}_spikein.sam"
    
    echo "Spike-in alignment complete for: $sample"
}

# Process all samples
for r1 in "$FASTQ_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    align_spikein "$sample"
done

echo "All spike-in alignments complete!"
```

**Expected Output**:
- `{sample}_spikein.bam`: Reads aligned to spike-in controls
- `{sample}_unmapped.fastq.gz`: Unmapped reads (for genome alignment)
- `{sample}_spikein_alignment.log`: Alignment statistics

---

### Step 3: Align to Human Genome

**Script**: `03_align_genome.sh`

```bash
#!/bin/bash
# Align spike-in unmapped reads to human genome (hg38)

set -e

UNMAPPED_DIR="/data/processed_data/spikein_alignment"
GENOME_REF="/data/references/hg38/hg38"
OUTPUT_DIR="/data/processed_data/bam"
THREADS=8

mkdir -p "$OUTPUT_DIR"

align_genome() {
    local sample=$1
    local r1="${UNMAPPED_DIR}/${sample}_unmapped.fastq.1.gz"
    local r2="${UNMAPPED_DIR}/${sample}_unmapped.fastq.2.gz"
    
    echo "=== Aligning $sample to hg38 ==="
    
    # Align to human genome
    bowtie2 \
        -x "$GENOME_REF" \
        -1 "$r1" \
        -2 "$r2" \
        --threads "$THREADS" \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        --maxins 1000 \
        -S "${OUTPUT_DIR}/${sample}_genome.sam" \
        2> "${OUTPUT_DIR}/${sample}_genome_alignment.log"
    
    # Convert to BAM, filter, sort
    samtools view -@ "$THREADS" -bS -q 10 "${OUTPUT_DIR}/${sample}_genome.sam" | \
        samtools sort -@ "$THREADS" -o "${OUTPUT_DIR}/${sample}.bam"
    
    # Index BAM file
    samtools index "${OUTPUT_DIR}/${sample}.bam"
    
    # Generate alignment statistics
    samtools flagstat "${OUTPUT_DIR}/${sample}.bam" > \
        "${OUTPUT_DIR}/${sample}_flagstat.txt"
    
    # Remove SAM file
    rm "${OUTPUT_DIR}/${sample}_genome.sam"
    
    echo "Genome alignment complete for: $sample"
}

# Process all samples
for unmapped in "$UNMAPPED_DIR"/*_unmapped.fastq.1.gz; do
    sample=$(basename "$unmapped" _unmapped.fastq.1.gz)
    align_genome "$sample"
done

echo "All genome alignments complete!"
```

**Expected Output**:
- `{sample}.bam`: Aligned reads (MAPQ ≥ 10)
- `{sample}.bam.bai`: BAM index
- `{sample}_flagstat.txt`: Alignment statistics

---

### Step 4: Calculate Spike-In Normalization Factors

**Script**: `04_calculate_normalization.R`

```r
#!/usr/bin/env Rscript
# Calculate spike-in normalization factors
# Based on Wilson et al. 2022 methodology

library(GenomicRanges)
library(Rsamtools)
library(tidyverse)

# Directories
spikein_dir <- "/data/processed_data/spikein_alignment"
output_dir <- "/data/results/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get list of samples
bam_files <- list.files(spikein_dir, pattern = "_spikein.bam$", full.names = TRUE)
samples <- gsub("_spikein.bam", "", basename(bam_files))

cat("Processing", length(samples), "samples\n\n")

# Function to count spike-in aligned reads
count_spikein_reads <- function(bam_file) {
  param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
  aln <- scanBam(bam_file, param = param)
  return(length(aln[[1]]$qname))
}

# Count reads for each sample
spikein_counts <- tibble(
  sample = samples,
  spikein_reads = sapply(bam_files, count_spikein_reads)
)

# Calculate normalization factors
# Method 1: Reads per million (RPM) normalization
spikein_counts <- spikein_counts %>%
  mutate(
    total_reads = spikein_reads,  # In practice, add genome-aligned reads
    rpm_factor = spikein_reads / 1e6,
    # Spike-in normalization factor (relative to median)
    median_spikein = median(spikein_reads),
    norm_factor = spikein_reads / median_spikein,
    # Scaling factor for downstream analysis
    scale_factor = 1 / norm_factor
  )

# Print summary
cat("\nSpike-in alignment summary:\n")
print(spikein_counts)

# Save normalization factors
write_csv(spikein_counts, file.path(output_dir, "spikein_normalization_factors.csv"))

# Visualization: Spike-in read counts
library(ggplot2)

p1 <- ggplot(spikein_counts, aes(x = reorder(sample, spikein_reads), 
                                  y = spikein_reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = median(spikein_counts$spikein_reads), 
             linetype = "dashed", color = "red") +
  labs(
    title = "Spike-In Aligned Reads per Sample",
    x = "Sample",
    y = "Number of Reads",
    caption = "Red line = median"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "../figures/spikein_reads.pdf"), 
       p1, width = 10, height = 6)

# Visualization: Normalization factors
p2 <- ggplot(spikein_counts, aes(x = reorder(sample, norm_factor), 
                                  y = norm_factor)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Spike-In Normalization Factors",
    x = "Sample",
    y = "Normalization Factor (relative to median)",
    caption = "Values > 1 indicate higher spike-in recovery"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "../figures/normalization_factors.pdf"), 
       p2, width = 10, height = 6)

cat("\nNormalization factors saved to:", 
    file.path(output_dir, "spikein_normalization_factors.csv"), "\n")
cat("Plots saved to:", file.path(output_dir, "../figures/"), "\n")
```

**Expected Output**:
- `spikein_normalization_factors.csv`: Table with normalization factors
- `spikein_reads.pdf`: Bar plot of spike-in read counts
- `normalization_factors.pdf`: Bar plot of normalization factors

---

### Step 5: Generate Normalized Coverage Tracks

**Script**: `05_generate_coverage.sh`

```bash
#!/bin/bash
# Generate normalized coverage tracks (bedGraph format)

set -e

BAM_DIR="/data/processed_data/bam"
NORM_FILE="/data/results/tables/spikein_normalization_factors.csv"
OUTPUT_DIR="/data/processed_data/bedgraph"
THREADS=4

mkdir -p "$OUTPUT_DIR"

# Function to generate normalized coverage
generate_coverage() {
    local sample=$1
    local bam="${BAM_DIR}/${sample}.bam"
    local scale_factor=$2
    
    echo "Generating coverage for: $sample (scale factor: $scale_factor)"
    
    # Generate coverage with bedtools
    bedtools genomecov \
        -ibam "$bam" \
        -bg \
        -scale "$scale_factor" \
        > "${OUTPUT_DIR}/${sample}_raw.bedGraph"
    
    # Sort bedGraph
    sort -k1,1 -k2,2n "${OUTPUT_DIR}/${sample}_raw.bedGraph" > \
        "${OUTPUT_DIR}/${sample}.bedGraph"
    
    # Remove unsorted file
    rm "${OUTPUT_DIR}/${sample}_raw.bedGraph"
    
    # Convert to bigWig (optional, requires bedGraphToBigWig)
    if command -v bedGraphToBigWig &> /dev/null; then
        bedGraphToBigWig \
            "${OUTPUT_DIR}/${sample}.bedGraph" \
            /data/references/hg38/hg38.chrom.sizes \
            "${OUTPUT_DIR}/${sample}.bigWig"
    fi
    
    echo "Coverage generated for: $sample"
}

# Read normalization factors and generate coverage
tail -n +2 "$NORM_FILE" | while IFS=, read -r sample spikein_reads total_reads \
    rpm_factor median_spikein norm_factor scale_factor; do
    
    if [ -f "${BAM_DIR}/${sample}.bam" ]; then
        generate_coverage "$sample" "$scale_factor"
    fi
done

echo "All coverage tracks generated!"
```

**Expected Output**:
- `{sample}.bedGraph`: Normalized coverage track
- `{sample}.bigWig`: Normalized coverage in bigWig format (optional)

---

### Step 6: Window-Based Quantification

**Script**: `06_window_quantification.R`

```r
#!/usr/bin/env Rscript
# Window-based quantification of methylation signal
# Following cfMeDIP-seq standard approach

library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(tidyverse)

# Parameters
WINDOW_SIZE <- 300  # bp
STEP_SIZE <- 300    # Non-overlapping windows

# Directories
bam_dir <- "/data/processed_data/bam"
norm_file <- "/data/results/tables/spikein_normalization_factors.csv"
output_dir <- "/data/processed_data/counts"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load normalization factors
norm_factors <- read_csv(norm_file)

# Define genomic windows
# Using chromosome 1-22, X, Y
chrom_info <- read.table("/data/references/hg38/hg38.chrom.sizes", 
                         col.names = c("chr", "length"))
chrom_info <- chrom_info %>%
  filter(chr %in% paste0("chr", c(1:22, "X", "Y")))

# Create tiling windows
windows_list <- lapply(1:nrow(chrom_info), function(i) {
  chr <- chrom_info$chr[i]
  len <- chrom_info$length[i]
  
  starts <- seq(1, len - WINDOW_SIZE + 1, by = STEP_SIZE)
  ends <- starts + WINDOW_SIZE - 1
  ends[ends > len] <- len
  
  GRanges(
    seqnames = chr,
    ranges = IRanges(start = starts, end = ends)
  )
})

windows <- do.call(c, windows_list)
cat("Created", length(windows), "genomic windows\n")

# Function to count reads in windows
count_reads_in_windows <- function(bam_file, windows, scale_factor) {
  param <- ScanBamParam(
    what = c("rname", "pos", "qwidth"),
    flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE)
  )
  
  aln <- readGAlignments(bam_file, param = param)
  aln_gr <- granges(aln)
  
  # Count overlaps
  counts <- countOverlaps(windows, aln_gr)
  
  # Apply normalization
  normalized_counts <- counts * scale_factor
  
  return(normalized_counts)
}

# Process all samples
count_matrix <- matrix(nrow = length(windows), 
                      ncol = nrow(norm_factors))
colnames(count_matrix) <- norm_factors$sample

for (i in 1:nrow(norm_factors)) {
  sample <- norm_factors$sample[i]
  scale <- norm_factors$scale_factor[i]
  bam_file <- file.path(bam_dir, paste0(sample, ".bam"))
  
  cat("Processing:", sample, "\n")
  count_matrix[, i] <- count_reads_in_windows(bam_file, windows, scale)
}

# Create count data frame
count_df <- as.data.frame(count_matrix)
count_df$chr <- as.character(seqnames(windows))
count_df$start <- start(windows)
count_df$end <- end(windows)

# Reorder columns
count_df <- count_df %>%
  select(chr, start, end, everything())

# Save count matrix
write_csv(count_df, file.path(output_dir, "window_counts_normalized.csv"))

cat("\nWindow quantification complete!\n")
cat("Output saved to:", file.path(output_dir, "window_counts_normalized.csv"), "\n")
```

**Expected Output**:
- `window_counts_normalized.csv`: Normalized read counts in 300bp windows

---

### Step 7: Statistical Analysis with MEDIPS

**Script**: `07_medips_analysis.R`

```r
#!/usr/bin/env Rscript
# MEDIPS analysis for differential methylation
# Incorporating spike-in normalization

library(MEDIPS)
library(GenomicRanges)
library(edgeR)
library(tidyverse)

# Directories
bam_dir <- "/data/processed_data/bam"
norm_file <- "/data/results/tables/spikein_normalization_factors.csv"
output_dir <- "/data/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load normalization factors
norm_factors <- read_csv(norm_file)

# Sample metadata (modify based on your experimental design)
# Example: cancer vs. normal
sample_groups <- tibble(
  sample = norm_factors$sample,
  group = c(rep("cancer", 5), rep("normal", 5))  # Adjust as needed
)

# Create MEDIPS set for each sample
medips_sets <- list()

for (i in 1:nrow(norm_factors)) {
  sample <- norm_factors$sample[i]
  bam_file <- file.path(bam_dir, paste0(sample, ".bam"))
  
  cat("Creating MEDIPS set for:", sample, "\n")
  
  medips_sets[[sample]] <- MEDIPS.createSet(
    file = bam_file,
    BSgenome = "BSgenome.Hsapiens.UCSC.hg38",
    extend = 300,
    shift = 0,
    uniq = 1e-3,
    window_size = 300,
    chr.select = paste0("chr", c(1:22, "X"))
  )
}

# Combine MEDIPS sets
cat("\nCreating count matrix...\n")
count_matrix <- MEDIPS.createROIset(
  genomicCoordinates = medips_sets[[1]]@genome_CF,
  ROIs = medips_sets
)

# Apply spike-in normalization factors
for (i in 1:ncol(count_matrix)) {
  sample <- norm_factors$sample[i]
  scale <- norm_factors$scale_factor[i]
  count_matrix[, i] <- count_matrix[, i] * scale
}

# Differential methylation analysis with edgeR
cat("\nPerforming differential methylation analysis...\n")

# Create DGEList
groups <- factor(sample_groups$group)
dge <- DGEList(counts = count_matrix, group = groups)

# Filter low-count windows
keep <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize (TMM)
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ groups)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit GLM
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)

# Extract results
dmr_results <- topTags(qlf, n = Inf)$table
dmr_results$chr <- seqnames(medips_sets[[1]]@genome_CF)[rownames(dmr_results)]
dmr_results$start <- start(medips_sets[[1]]@genome_CF)[rownames(dmr_results)]
dmr_results$end <- end(medips_sets[[1]]@genome_CF)[rownames(dmr_results)]

# Save results
write_csv(dmr_results, file.path(output_dir, "differential_methylation_results.csv"))

# Volcano plot
library(ggplot2)

p <- ggplot(dmr_results, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("gray", "red")) +
  labs(
    title = "Differential Methylation Analysis",
    x = "Log2 Fold Change",
    y = "-log10(FDR)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "figures/volcano_plot.pdf"), 
       p, width = 8, height = 6)

cat("\nAnalysis complete!\n")
cat("Results saved to:", file.path(output_dir, "differential_methylation_results.csv"), "\n")
```

**Expected Output**:
- `differential_methylation_results.csv`: DMR analysis results
- `volcano_plot.pdf`: Visualization of differential methylation

---

## Expected Outputs

### Directory Structure After Pipeline

```
/data/
├── processed_data/
│   ├── qc/
│   │   └── fastqc/              # Quality control reports
│   ├── spikein_alignment/
│   │   ├── *_spikein.bam        # Spike-in aligned reads
│   │   └── *_unmapped.fastq.gz  # Unmapped reads for genome alignment
│   ├── bam/
│   │   ├── *.bam                # Genome-aligned reads
│   │   └── *_flagstat.txt       # Alignment statistics
│   ├── bedgraph/
│   │   ├── *.bedGraph           # Normalized coverage tracks
│   │   └── *.bigWig             # BigWig format (optional)
│   └── counts/
│       └── window_counts_normalized.csv
│
└── results/
    ├── tables/
    │   ├── spikein_normalization_factors.csv
    │   └── differential_methylation_results.csv
    └── figures/
        ├── spikein_reads.pdf
        ├── normalization_factors.pdf
        └── volcano_plot.pdf
```

---

## Quality Control Metrics

### Key Metrics to Monitor

1. **FastQC Metrics**
   - Per base sequence quality > Q30
   - Adapter contamination < 5%
   - Sequence duplication levels

2. **Spike-In Alignment**
   - Spike-in alignment rate: 0.1-5% of total reads
   - Consistent across samples (CV < 30%)

3. **Genome Alignment**
   - Overall alignment rate > 80%
   - Properly paired reads > 95%
   - Mapping quality ≥ 10

4. **Library Complexity**
   - Number of unique reads
   - PCR duplication rate < 30%

5. **Coverage Distribution**
   - Uniform coverage across chromosomes
   - Enrichment in CpG islands

### Expected Spike-In Performance

Based on Wilson et al. 2022:
- Spike-in reads should be 0.5-2% of total
- CV of spike-in reads across samples < 20% indicates good batch consistency
- Normalization factors should cluster around 1.0 (±0.3)

---

## Running the Complete Pipeline

### Sequential Execution

```bash
# Navigate to project directory
cd ~/projects/wilson_spikein

# Copy scripts to project
cp /path/to/scripts/*.sh .
cp /path/to/scripts/*.R .

# Make scripts executable
chmod +x *.sh

# Run pipeline step by step
./01_fastqc.sh
./02_align_spikein.sh
./03_align_genome.sh
Rscript 04_calculate_normalization.R
./05_generate_coverage.sh
Rscript 06_window_quantification.R
Rscript 07_medips_analysis.R
```

### Parallel Processing

For large datasets, use GNU Parallel:

```bash
# Create sample list
ls /data/raw_data/fastq/*_R1.fastq.gz | \
    xargs -n 1 basename | \
    sed 's/_R1.fastq.gz//' > samples.txt

# Parallel alignment to spike-ins
parallel -j 4 --joblog spikein_align.log \
    './process_single_sample_spikein.sh {}' :::: samples.txt

# Parallel genome alignment
parallel -j 4 --joblog genome_align.log \
    './process_single_sample_genome.sh {}' :::: samples.txt
```

---

## Troubleshooting

### Common Issues

1. **Low spike-in alignment rate**
   - Check spike-in sequences are correct
   - Verify spike-in was added to samples

2. **High variation in normalization factors**
   - Indicates batch effects
   - Review IP protocol consistency
   - Check antibody lot numbers

3. **Memory errors during R processing**
   - Process chromosomes separately
   - Use `data.table` for large datasets
   - Increase instance RAM (upgrade to r5.2xlarge)

4. **Slow Bowtie2 alignment**
   - Use `--threads` parameter
   - Consider spot instances for batch processing

---

## References

Wilson GA, et al. (2022). Synthetic spike-in controls improve detection of enrichment-based chromatin modification assays. *Cell Reports Methods*, 2(12):100366.

Shen SY, et al. (2019). Preparation of cfMeDIP-seq libraries for methylome profiling of plasma cell-free DNA. *Nature Protocols*, 14(10):2749-2780.

---

**Pipeline Version**: 1.0  
**Last Updated**: January 2026  
**Compatibility**: cfMeDIP-seq data with spike-in controls
