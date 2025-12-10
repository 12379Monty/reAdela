# Differential Methylation Analysis Methods: A Comprehensive Technical Report

---

## Executive Summary

This report provides a comprehensive overview of computational methods for differential methylation analysis across three major DNA methylation profiling platforms: enrichment-based sequencing (MeDIP-seq/cfMeDIP-seq), bisulfite sequencing (WGBS/RRBS), and methylation arrays (Illumina 450K/EPIC). The report is organized into three parts: (1) a detailed catalog of analytical methods with publication highlights, (2) an examination of the input data structures that distinguish each platform, and (3) a review of benchmark studies that have evaluated method performance head-to-head.

---

# PART I: Differential Methylation Analysis Methods

## 1. Enrichment-Based Methods (MeDIP-seq / cfMeDIP-seq)

Enrichment-based methods use antibody-mediated immunoprecipitation to capture methylated DNA fragments. These methods require specialized algorithms to convert enrichment signals into absolute methylation estimates and to identify differentially methylated regions.

### 1.1 Batman (Bayesian Tool for Methylation Analysis)

**Publication:** Down TA et al. A Bayesian deconvolution strategy for immunoprecipitation-based DNA methylome analysis. *Nature Biotechnology* 2008;26(7):779-785.

**Key Highlights:**

- First cross-platform algorithm developed to estimate absolute DNA methylation levels from MeDIP profiles
- Uses a Bayesian deconvolution strategy to account for CpG density effects on immunoprecipitation enrichment
- Core principle: models the effect that varying CpG density has on MeDIP enrichment of DNA fragments
- The algorithm assumes the MeDIP signal scales linearly with the number of methylated CpG sites in a sequence
- Can be applied to both MeDIP-chip (microarray) and MeDIP-seq (sequencing) data
- Generated the first high-resolution whole-genome DNA methylome of any mammalian genome (human spermatozoa)

**Statistical Model:** Batman uses a statistical model where the complete MeDIP dataset can be represented as a probability distribution incorporating CpG coupling factors between probes and CpG dinucleotides. Standard Bayesian techniques using nested sampling are employed to infer the distribution of likely methylation states.

**Limitations:** Computationally intensive (can take several days to weeks per chromosome); requires consideration of copy number variation.

---

### 1.2 MEDIPS (Methylated DNA Immunoprecipitation Sequencing)

**Publication:** Lienhard M et al. MEDIPS: genome-wide differential coverage analysis of sequencing data derived from DNA enrichment experiments. *Bioinformatics* 2014;30(2):284-286.

**Key Highlights:**

- R/Bioconductor package for genome-wide MeDIP-seq analysis
- Implements linear regression to model the relationship between read density and methylation level
- Provides MeDIP-seq-specific quality control metrics including:
  - Saturation analysis to assess sequencing depth adequacy
  - CpG enrichment factor calculation to validate immunoprecipitation success
  - Reproducibility assessment across technical replicates
- Includes methodology for identifying differentially methylated regions (DMRs) between samples

**CpG Density Correction:** MEDIPS uses a coupling factor score to normalize for CpG density bias, though this linear model may not fully capture the saturation behavior observed at high CpG density regions.

**Limitation:** The original DMR calling algorithm requires an input (non-immunoprecipitated) sample in addition to the IP sample, effectively doubling sequencing costs.

---

### 1.3 QSEA (Quantitative Sequence Enrichment Analysis)

**Publication:** Lienhard M et al. QSEA—modelling of genome-wide DNA methylation from sequencing enrichment experiments. *Nucleic Acids Research* 2017;45(6):e44.

**Key Highlights:**

- Comprehensive workflow for modeling and quantifying MeDIP-seq data
- Introduces a Bayesian statistical model that transforms enrichment read counts to absolute methylation levels
- Uses a sigmoidal model (rather than linear) with empirical calibration to better capture saturation effects at high CpG density
- Enhances interpretability by enabling direct comparison with other methylation assays (e.g., bisulfite sequencing)
- Provides calibration strategies using either additional validation data or general assumptions about methylation patterns
- Available as a Bioconductor package

**Validation:** Successfully validated on a clinically relevant benchmark dataset from non-small cell lung cancer patients; retrieved well-known lung tumor methylation markers causative for gene expression changes.

**Key Innovation:** Sigmoidal model accounts for antibody binding saturation in CpG-rich regions; can work without bisulfite sequencing calibration data by leveraging assumptions about CpG island methylation patterns.

---

### 1.4 MeDEStrand

**Publication:** Xu J et al. MeDEStrand: an improved method to infer genome-wide absolute methylation levels from DNA enrichment data. *BMC Bioinformatics* 2018;19:540.

**Key Highlights:**

- Method for inferring genome-wide absolute methylation levels from DNA enrichment data
- Key innovation: processes reads for positive and negative DNA strands separately
- Uses a logistic regression (sigmoid) function to estimate and correct CpG bias
- Addresses asymmetric CpG methylation across DNA strands
- Demonstrated superior performance at high resolution (25, 50, and 100 bp bins)

**Validation:** Compared against MEDIPS, BayMeth, and QSEA on four independent datasets (GM12878, K562, foreskin fibroblasts, mammary epithelial cells) using RRBS as ground truth.

**Improvements over MEDIPS:**
1. Uses logistic regression model (sigmoid function) instead of linear regression for CpG density correction
2. Separately estimates and corrects CpG bias for positive and negative DNA strands

**Availability:** R package at https://github.com/jxu1234/MeDEStrand.git

---

### 1.5 MEDIPIPE

**Publication:** Zeng Y et al. MEDIPIPE: an automated and comprehensive pipeline for cfMeDIP-seq data quality control and analysis. *Bioinformatics* 2023;39(7):btad423.

**Key Highlights:**

- Automated end-to-end pipeline specifically designed for cfMeDIP-seq data
- Provides comprehensive quality control framework for cell-free methylated DNA analysis
- Developed using Snakemake for reproducible, containerized execution
- Supports various experimental settings through a single configuration file:
  - Single-end or paired-end sequencing
  - Spike-in controls
  - Unique molecular identifiers (UMIs)
- Integrates three methylation quantification methods:
  - MEDIPS (linear regression)
  - QSEA (sigmoidal model with empirical knowledge)
  - MEDStrand (stranded sigmoid model)

**Four Main Modules:**
1. Read trimming and QC
2. Processed read alignment and QC
3. DNA methylation quantification
4. Data aggregation and filtering

**Availability:** https://github.com/pughlab/MEDIPIPE

---

### 1.6 GLMnet-Based Classification (cfMeDIP-seq Cancer Detection)

**Publication:** Shen SY et al. Sensitive tumour detection and classification using plasma cell-free DNA methylomes. *Nature* 2018;563(7732):579-583.

**Key Highlights:**

- Foundational approach for cfMeDIP-seq-based cancer classification
- Workflow: cfMeDIP-seq paired-end data reduced to 300 bp genomic windows mapping to CpG islands, shores, shelves, and FANTOM5 enhancers
- Feature selection: Top 1,000 most variable fragments between cancer patients and cancer-free controls
- Classification using GLMnet (elastic net regularized generalized linear model)
- Also explored Random Forest as alternative classifier
- Demonstrated ability to discriminate between different cancer types based on plasma cfDNA methylation

**Data Availability:** Machine learning models and intermediate data objects publicly available on Zenodo.

---

### 1.7 Alternative Classification Approaches for cfMeDIP-seq

**Publication:** Halla-aho V & Lähdesmäki H. Probabilistic modeling methods for cell-free DNA methylation based cancer classification. *BMC Bioinformatics* 2022;23:138.

**Key Highlights:**

- Comprehensive comparison of statistical methods for cfMeDIP-seq cancer classification
- Evaluated multiple feature selection methods: t-tests, moderated t-tests, Fisher's exact test
- Evaluated feature generation approaches: PCA, Iterative Supervised PCA (ISPCA)
- Evaluated classification methods: GLMnet, logistic regression with sparsity-promoting priors

**Key Findings:**
- Many methods perform as well as or better than GLMnet in certain contexts
- Performance varies across sequencing depths, cancer types, and study cohorts
- Fisher's exact test and ISPCA showed robust feature selection
- Simple logistic regression using counts of hyper/hypo-methylated regions performed well

**Recommended Robust Methods:**
- Feature selection: Fisher's exact test, ISPCA
- Classification: Logistic regression with hyper/hypo-methylated region counts as features

---

### 1.8 MethRaFo

**Publication:** Ding J & Bar-Joseph Z. MethRaFo: MeDIP-seq methylation estimate using a Random Forest Regressor. *Bioinformatics* 2017;33(21):3477-3479.

**Key Highlights:**

- Uses Random Forest regression for MeDIP-seq methylation correction
- Trained on datasets profiled with both BS-seq and MeDIP-seq
- Achieved ~4-fold decrease in runtime compared to MEDIPS and BayMeth
- Increased accuracy by up to 20% over prior methods
- Uses nucleotide-level RPKM and local CpG density as features

**Performance Comparison:** Batman could not complete within 1 week on test data; MEDIPS and BayMeth took ~1 hour each with ~12GB RAM usage.

---

## 2. Bisulfite Sequencing Methods (WGBS / RRBS)

Bisulfite sequencing provides single-nucleotide resolution methylation data by converting unmethylated cytosines to uracil. These methods require statistical approaches that account for the binomial/beta-binomial nature of the count data.

### 2.1 DSS (Dispersion Shrinkage for Sequencing)

**Publication:** Feng H et al. A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. *Nucleic Acids Research* 2014;42(8):e69.

**Key Highlights:**

- Bioconductor package implementing Bayesian hierarchical model for DML detection
- Uses beta-binomial distribution to model BS-seq count data
- Employs empirical Bayes approach for dispersion shrinkage across all CpG sites
- Implements rigorous Wald test for hypothesis testing
- Handles both DML (differentially methylated locus) and DMR detection

**Extension - DSS-general (2016):**
- Supports complex experimental designs with multiple groups and covariates
- Uses arcsine link function for computational efficiency (50x faster than canonical GLM)

**Statistical Innovation:** Prior distribution constructed from whole genome; arcsine link function breaks variance-mean dependence, allowing weighted least squares instead of iterative GLM fitting.

---

### 2.2 methylKit

**Publication:** Akalin A et al. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. *Genome Biology* 2012;13:R87.

**Key Highlights:**

- Comprehensive R package for genome-wide DNA methylation profile analysis
- Provides functions for data import, quality control, clustering, visualization, and differential methylation analysis
- Statistical approaches:
  - Fisher's exact test (for comparisons without replicates)
  - Logistic regression with overdispersion correction (when replicates available)
- Automatically selects appropriate statistical test based on sample size
- P-values corrected using Benjamini-Hochberg FDR method
- Supports regional analysis through user-defined genomic windows or annotations

**Usage Guidelines:**
- Single sample per group → Fisher's exact test
- Multiple samples per group → Logistic regression
- Demonstrated good performance at low sequencing depths

---

### 2.3 BSmooth

**Publication:** Hansen KD et al. BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions. *Genome Biology* 2012;13:R83.

**Key Highlights:**

- Pipeline designed specifically for whole genome bisulfite sequencing (WGBS) data
- Key innovation: local likelihood smoothing of methylation estimates
- Uses smoothed methylation values across genomic regions
- Subsequent differential methylation testing via linear regression on smoothed values
- Performs analysis on methylation profiles rather than individual CpG sites
- Can identify large differentially methylated regions

**Considerations:**
- Smoothing may reduce sensitivity for detecting DMCs at individual sites
- Particularly suited for identifying broad regional methylation differences
- Retains sign of differential methylation during smoothing, which may cause signal cancellation at direction change boundaries
- Requires at least 3 replicates per condition

---

### 2.4 BiSeq

**Publication:** Hebestreit K et al. Detection of significantly differentially methylated regions in targeted bisulfite sequencing data. *Bioinformatics* 2013;29(13):1647-1653.

**Key Highlights:**

- Designed for detection of significantly differentially methylated regions in targeted bisulfite sequencing
- Uses beta regression framework on smoothed methylation levels
- Groups nearby CpG sites into clusters for regional analysis
- Performs differential methylation testing on clustered regions
- Showed highest sensitivity on data with small sequencing coverage due to smoothing approach

**Caution:** In comprehensive evaluations, BiSeq showed high false positive rates in null model testing, suggesting potential for overdetection.

---

### 2.5 RADMeth

**Publication:** Dolzhenko E & Smith AD. A beta-binomial mixture model for molecular count data. *Bioinformatics* 2014.

**Key Highlights:**

- Regression analysis of differential methylation
- Uses beta-binomial regression framework
- Designed to handle multiple covariates in experimental design
- Demonstrated slightly better performance than methylKit and DSS at low sequencing depth
- Showed higher sensitivity overall in DMC detection
- Achieves good balance between sensitivity and specificity

---

### 2.6 MOABS (Model-Based Analysis of Bisulfite Sequencing)

**Publication:** Sun D et al. MOABS: model based analysis of bisulfite sequencing data. *Genome Biology* 2014;15:R38.

**Key Highlights:**

- Comprehensive suite for BS-seq data analysis
- Uses beta-binomial assumption with bimodal prior distribution
- Introduces "credible methylation difference" metric combining:
  - Biological significance (magnitude of methylation difference)
  - Statistical significance (confidence in the difference)
- Employs empirical Bayes approach for parameter estimation
- When biological replicates available, uses maximum likelihood estimation
- Provides single unified metric for ranking DMRs

---

### 2.7 metilene

**Publication:** Jühling F et al. metilene: fast and sensitive calling of differentially methylated regions from bisulfite sequencing data. *Genome Research* 2016;26(2):256-262.

**Key Highlights:**

- Fast and sensitive DMR caller
- Uses a binary segmentation algorithm combined with a two-dimensional statistical test
- Does not require pre-defined genomic regions
- Merges CpGs based on actual genomic position
- Not effective for single-CpG analysis but performs well for regional analysis
- More accurately detects DMR boundaries compared to window-based methods

---

### 2.8 HMM-Based Methods (HMM-DM and HMM-Fisher)

**Publication:** Yu X & Sun S. HMM-DM: identifying differentially methylated regions using a hidden Markov model. *Statistical Applications in Genetics and Molecular Biology* 2016;15(1):69-81.

**Key Highlights:**

- Uses Hidden Markov Models along the genome to infer DMRs
- HMM-Fisher combines adjacent CpG sites (within 100 bases) during Fisher's exact test
- Relatively higher sensitivity and lower false positive rates than other methods
- Particularly effective for DMRs with large variation
- Suitable for small sample sizes

---

## 3. Array-Based Methods (450K / EPIC)

Illumina methylation arrays measure methylation at predefined CpG sites. Analysis methods must account for the beta-distributed nature of methylation values and probe-level technical considerations.

### 3.1 limma (Linear Models for Microarray Data)

**Publication:** Smyth GK. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology* 2004;3:Article3.

**Key Highlights:**

- Gold-standard statistical framework for differential methylation analysis
- Uses empirical Bayes framework based on Gaussian model theory
- Computes moderated t-statistics by borrowing information across probes
- Supports complex experimental designs through design matrices

**Standard Workflow for 450K/EPIC Arrays:**
1. Convert beta values to M-values (logit transformation)
2. Fit linear models using design matrix
3. Apply empirical Bayes shrinkage
4. Adjust for multiple testing

**Critical Consideration:** M-values (logit-transformed beta values) are recommended for statistical analysis due to better statistical properties, while beta values are more interpretable for biological effect sizes.

---

### 3.2 DMRcate

**Publication:** Peters TJ et al. De novo identification of differentially methylated regions in the human genome. *Epigenetics & Chromatin* 2015;8:6.

**Key Highlights:**

- Novel method for de novo identification of differentially methylated regions
- Uses tunable kernel (Gaussian) smoothing of differential methylation signal
- Key features:
  - Agnostic to genomic annotation (doesn't require predefined regions)
  - Agnostic to direction of methylation change within regions
  - Removes bias from irregularly spaced methylation sites
- Uses unsigned weights (squared t-statistics) to avoid signal cancellation
- Assigns significance to DMRs via comparison to null model
- Fast computational performance (no permutations required)

**Performance:** Outperformed Bumphunter and Probe Lasso in predictive performance; commensurate performance with comb-p.

---

### 3.3 Bumphunter

**Publication:** Jaffe AE et al. Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies. *International Journal of Epidemiology* 2012;41(1):200-209.

**Key Highlights:**

- Identifies differentially methylated regions by "bumphunting" approach
- Smooths methylation values across genomic regions
- Uses permutation testing for significance assignment
- Implemented in minfi Bioconductor package
- Retains sign of differential methylation during smoothing

**Considerations:**
- Computationally intensive due to permutation-based significance testing
- May have reduced sensitivity at boundaries where methylation direction changes abruptly
- Works well with high-density probe coverage

---

### 3.4 comb-p

**Publication:** Pedersen BS et al. Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values. *Bioinformatics* 2012;28(22):2986-2988.

**Key Highlights:**

- Combines spatially correlated p-values to identify regions of enrichment
- Uses Stouffer-Liptak correction for spatial autocorrelation
- Can be applied to any probe-wise analysis results
- Does not require specific input data format

**Note:** Has shown high false positive rates in some benchmark studies.

---

### 3.5 missMethyl (GOmeth and GOregion)

**Publication:** Phipson B et al. missMethyl: an R package for analyzing data from Illumina's HumanMethylation450 platform. *Bioinformatics* 2016;32(2):286-288.

**Key Highlights:**

- Bioconductor package extending limma for methylation array analysis
- Key contributions:
  - RUV (Remove Unwanted Variation) integration for batch effect correction
  - GOmeth: unbiased gene set testing for differentially methylated CpGs
  - GOregion: gene set testing for differentially methylated regions

**GOmeth Innovation:** Addresses bias in gene set testing by accounting for varying numbers of probes per gene (genes with more probes are more likely to appear significant by chance).

---

### 3.6 coMethDMR

**Publication:** Gomez L et al. coMethDMR: accurate identification of co-methylated and differentially methylated regions in epigenome-wide association studies with continuous phenotypes. *Nucleic Acids Research* 2019;47(17):e98.

**Key Highlights:**

- Identifies co-methylated regions first, then tests for differential methylation
- Uses correlation-based clustering to define regions
- Particularly well-suited for continuous phenotypes
- Showed well-controlled false positive rates (~5%) in benchmark studies

---

# PART II: Input Data Structures

The form of input data to differential methylation analysis modules differs fundamentally between platform types, which dictates the appropriate statistical models for each.

## 1. Enrichment-Based Methods (MeDIP-seq / cfMeDIP-seq)

### Input Data Structure: Binned Read Counts or Read Density

The input to differential methylation analysis is **read counts aggregated into genomic bins/windows**, not individual CpG-level measurements. This is because enrichment-based methods capture fragments containing methylated CpGs rather than measuring methylation at specific positions.

**Typical Data Format:**

| Column | Description |
|--------|-------------|
| Chromosome | Genomic chromosome |
| Start | Bin start position |
| End | Bin end position |
| Read Count (or RPKM) | Number of reads mapping to this bin |

**Typical bin sizes:** 100 bp, 300 bp, or 500 bp windows

**Example Data:**
```
chr1    10000    10100    45      # 45 reads in this 100bp bin
chr1    10100    10200    12      # 12 reads in next bin
chr1    10200    10300    89      # 89 reads (perhaps CpG island)
```

**Key Characteristics:**
- Read density is proportional to methylation level, but confounded by CpG density
- Higher CpG density regions attract more antibody binding → more reads (even at same methylation %)
- This is why CpG density correction (via Batman, MEDIPS, QSEA, etc.) is essential before differential analysis

**Appropriate Statistical Models:**
- Linear regression (MEDIPS)
- Negative binomial models (similar to RNA-seq)
- Poisson-based models
- Normalized enrichment scores after CpG density correction

---

## 2. Bisulfite Sequencing Methods (WGBS / RRBS)

### Input Data Structure: Per-CpG Count Data (Binomial Observations)

The input is **count data at each individual CpG site**: the number of reads showing methylation (unconverted C) versus the total number of reads covering that position.

**Typical Data Format:**

| Column | Description |
|--------|-------------|
| Chromosome | Genomic chromosome |
| Position | Genomic coordinate of the CpG |
| N (Total) | Total reads covering this CpG |
| X (Methylated) | Reads showing methylation (C retained) |

**Example Data (DSS input format):**
```
chr     pos         N       X
chr1    10542       18      15
chr1    10563       24      21
chr1    10571       31      28
chr1    10590       17      14
```

**Key Characteristics:**
- Each CpG is an independent binomial observation: X ~ Binomial(N, p) where p = true methylation proportion
- Methylation ratio = X/N (but using raw counts preserves uncertainty information)
- Coverage (N) varies across sites and samples
- Low coverage sites have high variance in methylation estimates

**Appropriate Statistical Models:**
- Beta-binomial distribution (DSS, MOABS) - accounts for biological overdispersion
- Binomial/logistic regression (methylKit with replicates)
- Fisher's exact test (methylKit without replicates) - comparing 2×2 tables of methylated/unmethylated counts

---

## 3. Array-Based Methods (450K / EPIC)

### Input Data Structure: Beta Values or M-Values

The input is **continuous methylation measurements** at each probe, typically expressed as beta values (proportion methylated) or M-values (logit-transformed beta values).

**Beta Value:** β = Methylated Signal / (Methylated Signal + Unmethylated Signal + offset)

**M-Value:** M = log2(β / (1 - β))

**Typical Data Format:**

| Probe ID | Sample1 | Sample2 | Sample3 | ... |
|----------|---------|---------|---------|-----|
| cg00000029 | 0.85 | 0.82 | 0.88 | ... |
| cg00000108 | 0.12 | 0.15 | 0.11 | ... |
| cg00000165 | 0.45 | 0.48 | 0.42 | ... |

**Key Characteristics:**
- Beta values bounded between 0 and 1
- M-values unbounded, approximately normally distributed
- Fixed probe locations (no coverage variation)
- ~450,000 probes (450K) or ~850,000 probes (EPIC)

**Appropriate Statistical Models:**
- Linear models with empirical Bayes (limma) on M-values
- t-tests on M-values or beta values
- Moderated statistics with variance shrinkage

---

## Summary Comparison of Input Data Structures

| Aspect | Enrichment-Based | Bisulfite Sequencing | Arrays |
|--------|------------------|---------------------|--------|
| **Unit of measurement** | Genomic bin (100-500 bp) | Individual CpG site | Individual probe |
| **Data type** | Read counts/density (continuous-like) | Binomial counts (discrete) | Continuous (beta/M-values) |
| **Raw input** | Reads per bin | (Methylated reads, Total reads) per CpG | Intensity ratios |
| **Resolution** | ~100-500 bp | Single nucleotide | Single probe |
| **Key confounder** | CpG density | Sequencing coverage | Probe bias (Type I vs II) |
| **Appropriate models** | Negative binomial, Poisson, linear regression | Beta-binomial, binomial, logistic regression | Linear models, empirical Bayes |
| **Variance structure** | Overdispersion from technical/biological variation | Overdispersion beyond binomial expectation | Heteroscedasticity |

**Critical Implication:** Methods designed for one platform type cannot be directly applied to another due to fundamentally different data generating processes and likelihood functions.

---

# PART III: Benchmark Studies and Method Comparisons

## 1. Enrichment-Based Methods (MeDIP-seq)

There is **no single comprehensive benchmark paper** comparing all MeDIP-seq differential methylation methods head-to-head. However, several papers provide partial comparisons:

### 1.1 QSEA Benchmark (Lienhard et al. 2017)

**Publication:** *Nucleic Acids Research* 45:e44

**Methods Compared:** QSEA vs BayMeth vs MeSiC

**Validation Approach:** Comparison against bisulfite sequencing on both in vitro samples and in vivo samples from non-small cell lung cancer (NSCLC) patients

**Key Findings:**
- QSEA outperformed both BayMeth and MeSiC
- QSEA particularly suited for situations where no additional calibration experiments are available
- Successfully retrieved well-known lung tumor methylation markers

---

### 1.2 MeDEStrand Benchmark (Xu et al. 2018)

**Publication:** *BMC Bioinformatics* 19:540

**Methods Compared:** MeDEStrand vs MEDIPS vs BayMeth vs QSEA

**Validation Approach:** Comparison against RRBS data on 4 cell types (GM12878, K562, foreskin fibroblasts, mammary epithelial cells)

**Key Findings:**
- MeDEStrand showed best performance at high resolution (25, 50, 100 bp bins)
- Strand-specific processing improved accuracy
- Sigmoid function better captured saturation effects than linear models

---

### 1.3 MethRaFo Benchmark (Ding & Bar-Joseph 2017)

**Publication:** *Bioinformatics* 33:3477-3479

**Methods Compared:** MethRaFo (Random Forest) vs MEDIPS vs BayMeth vs Batman

**Validation Approach:** Comparison against BS-seq on Roadmap Epigenomics data (brain, cortex, penis foreskin samples)

**Key Findings:**
- MethRaFo achieved ~4x faster runtime with ~20% accuracy improvement
- Batman could not complete within 1 week on genome-wide data
- MEDIPS and BayMeth took ~1 hour with ~12GB RAM usage
- Random Forest approach provided efficient alternative to Bayesian methods

---

## 2. Bisulfite Sequencing Methods (WGBS/RRBS)

The bisulfite sequencing field has the most mature benchmarking literature with several comprehensive comparison studies.

### 2.1 Comprehensive Evaluation (2021)

**Publication:** Park Y et al. Comprehensive Evaluation of Differential Methylation Analysis Methods for Bisulfite Sequencing Data. *International Journal of Environmental Research and Public Health* 2021;18:7975.

**Methods Compared:** methylKit, DSS, BSmooth, BiSeq, RADMeth, methylSig, Fisher's exact test, metilene

**Evaluation Metrics:**
- True positive rate at various sequencing depths (5x, 10x, 15x, 20x, 25x)
- Performance with varying replicate numbers (1-5 per condition)
- False positive rate under null model
- DMR boundary accuracy

**Key Findings:**

| Finding | Details |
|---------|---------|
| Overall | No single method consistently ranked first across all benchmarks |
| Low depth | RADMeth performed slightly better at low sequencing depth |
| Few replicates | DSS and methylKit had higher true positive rates with few replicates |
| Smoothing | Smoothing did NOT improve DMC detection accuracy, even for low-depth data |
| False positives | BiSeq showed very high FP rates (15,498 vs 0-37 for others in null model) |
| Threshold | methylKit, DSS, RADMeth performed well when coverage ≥15x or replicates ≥3 |
| DMR boundaries | metilene more accurately detected DMR boundaries than window-based methods |

**Null Model Results (False Positives):**
- methylSig: 0
- Fisher's exact test: 0
- DSS: 4
- methylKit: 37
- RADMeth: 36
- BiSeq: 15,498

---

### 2.2 Five-Method Comparison (Yu & Sun 2016)

**Publication:** Yu X & Sun S. Comparing five statistical methods of differential methylation identification using bisulfite sequencing data. *Statistical Applications in Genetics and Molecular Biology* 2016;15(1):69-81.

**Methods Compared:** methylKit, BSmooth, BiSeq, HMM-DM, HMM-Fisher

**Key Findings:**
- Parameter settings largely affect accuracy; modified settings often yield better results than defaults
- All methods perform better on long DMRs with small within-group variation
- Low concordance between methods due to different underlying approaches
- HMM-DM and HMM-Fisher showed higher sensitivity and lower false positive rates
- BiSeq best preserved raw methylation signals among smoothing methods

**Recommendations:** Select methods based on data characteristics and the advantages of each method.

---

### 2.3 Comprehensive Survey (Shafi et al. 2018)

**Publication:** Shafi A et al. A survey of the approaches for identifying differential methylation using bisulfite sequencing data. *Briefings in Bioinformatics* 2018;19(5):737-753.

**Methods Reviewed:** 20+ methods

**Classification Schemes:**

**By Input Data Type:**
- Count data: methylKit, eDMR, DSS, DSS-single, DSS-general, MOABS, RADMeth, MethylSig, MACAU, GetisDMR, ComMet
- Ratio data: BSmooth, BiSeq, qDMR, CpG_MPs, SMART, HMM-Fisher, HMM-DM, COHCAP, metilene
- Both: DMAP, swDMR

**Recommendations by Scenario:**

| Scenario | Recommended Methods |
|----------|-------------------|
| Small sample size | DSS, MethylSig, HMM-Fisher |
| Multiple covariates | methylKit, RADMeth, DSS-general, BiSeq, eDMR, MACAU, GetisDMR |
| Single sample patterns | QDMR, CpG_MPs, HMM-Fisher |
| Cell type-specific marks | SMART |

---

### 2.4 Strategies for Bisulfite Analysis (Akalin et al. 2017)

**Publication:** Akalin A et al. Strategies for analyzing bisulfite sequencing data. *bioRxiv* preprint.

**Methods Compared:** limma, DSS, BSmooth, methylKit (with/without overdispersion)

**Key Findings:**
- methylKit without overdispersion: highest F-score but lowest specificity (exploratory analysis)
- methylKit with overdispersion: higher specificity, second-highest F-score (balanced approach)
- limma: highest specificity but extremely few true positives
- BSmooth: not suitable for data without smooth methylation profiles

**Recommendation:** methylKit without overdispersion for exploratory analysis; methylKit with overdispersion, DSS, or limma when limiting false positives is critical.

---

## 3. Array-Based Methods (450K/EPIC)

### 3.1 False Positive Rate Evaluation (2022)

**Publication:** Gatev E et al. An evaluation of the genome-wide false positive rates of common methods for identifying differentially methylated regions using Illumina methylation arrays. *Epigenetics* 2022.

**Methods Compared:** DMRcate, Bumphunter, comb-p, mCSEA, coMethDMR

**Evaluation Approach:** Genome-wide null simulations with simulated phenotypes generated independently of methylation data; tested on both 450K and EPIC arrays with continuous and dichotomous phenotypes.

**Key Findings:**

| Method | False Positive Rate | Notes |
|--------|-------------------|-------|
| coMethDMR | ~5% (well-controlled) | Except for skewed continuous phenotypes |
| DMRcate | Generally well-controlled | Variable on EPIC; good on 450K |
| mCSEA | ≥0.096 | Elevated false positives |
| comb-p | >0.34 | High false positive rates |
| Bumphunter | 0.35-0.95 | Very high false positive rates |

---

### 3.2 DMRcate Benchmark (Peters et al. 2015)

**Publication:** Peters TJ et al. De novo identification of differentially methylated regions in the human genome. *Epigenetics & Chromatin* 2015;8:6.

**Methods Compared:** DMRcate vs Bumphunter vs Probe Lasso vs comb-p

**Evaluation Approach:** Both simulated and real 450K data; compared sensitivity, specificity, and runtime

**Key Findings:**
- DMRcate outperformed Bumphunter and Probe Lasso in predictive performance
- Commensurate performance with comb-p
- DMRcate much faster (no permutations required)
- Bumphunter may suffer signal cancellation at direction change boundaries due to signed smoothing

---

### 3.3 Microarray and NGS Quality Assessment (2023)

**Publication:** Abramov S et al. Assessing the Differential Methylation Analysis Quality for Microarray and NGS Platforms. *International Journal of Molecular Sciences* 2023;24:8591.

**Methods Compared:** 
- Arrays: limma, T-test, dmpFinder (with/without variance shrinkage)
- NGS: methylKit, DSS, BSmooth, BiSeq, RADMeth, HMM-DM

**Key Findings:**
- For arrays: methods produced largely overlapping signatures; most had well-controlled false positive rates
- For NGS: much lower concordance between methods; signatures varied dramatically
- DSS with smoothing and RADMeth had highest recall on simulated data
- methylKit produced shorter signatures than other methods

---

## 4. Summary of Benchmark Findings

### Key Benchmark Papers by Platform

| Platform | Paper | Year | Methods Compared | Key Recommendation |
|----------|-------|------|------------------|-------------------|
| MeDIP-seq | QSEA (Lienhard) | 2017 | QSEA, BayMeth, MeSiC | QSEA when no calibration data |
| MeDIP-seq | MeDEStrand (Xu) | 2018 | MeDEStrand, MEDIPS, QSEA, BayMeth | MeDEStrand for high resolution |
| MeDIP-seq | MethRaFo (Ding) | 2017 | MethRaFo, MEDIPS, BayMeth, Batman | MethRaFo for speed |
| BS-seq | Comprehensive Eval | 2021 | 8 methods | RADMeth/methylKit for low depth |
| BS-seq | Yu & Sun | 2016 | 5 methods | HMM-based methods for sensitivity |
| BS-seq | Survey (Shafi) | 2018 | 20+ methods | Context-dependent selection |
| Arrays | False Positive Eval | 2022 | 5 DMR methods | coMethDMR for controlled FPR |
| Arrays | DMRcate (Peters) | 2015 | 4 DMR methods | DMRcate for speed + accuracy |

### Cross-Platform Observations

1. **No universal best method:** Performance depends on data characteristics, sample size, and analysis goals
2. **Smoothing trade-offs:** Smoothing helps regional detection but may reduce single-CpG sensitivity
3. **Replicates matter more than depth:** For bisulfite sequencing, having more biological replicates is often more important than increased sequencing depth
4. **False positive control varies:** Some methods (e.g., BiSeq, Bumphunter) show elevated false positive rates in certain conditions
5. **Method concordance is low:** Different methods often identify different DMRs/DMCs, suggesting complementary use may be beneficial

---

# References

1. Down TA et al. A Bayesian deconvolution strategy for immunoprecipitation-based DNA methylome analysis. *Nat Biotechnol*. 2008;26(7):779-785.

2. Lienhard M et al. MEDIPS: genome-wide differential coverage analysis of sequencing data derived from DNA enrichment experiments. *Bioinformatics*. 2014;30(2):284-286.

3. Lienhard M et al. QSEA—modelling of genome-wide DNA methylation from sequencing enrichment experiments. *Nucleic Acids Res*. 2017;45(6):e44.

4. Xu J et al. MeDEStrand: an improved method to infer genome-wide absolute methylation levels from DNA enrichment data. *BMC Bioinformatics*. 2018;19:540.

5. Zeng Y et al. MEDIPIPE: an automated and comprehensive pipeline for cfMeDIP-seq data quality control and analysis. *Bioinformatics*. 2023;39(7):btad423.

6. Shen SY et al. Sensitive tumour detection and classification using plasma cell-free DNA methylomes. *Nature*. 2018;563(7732):579-583.

7. Halla-aho V & Lähdesmäki H. Probabilistic modeling methods for cell-free DNA methylation based cancer classification. *BMC Bioinformatics*. 2022;23:138.

8. Ding J & Bar-Joseph Z. MethRaFo: MeDIP-seq methylation estimate using a Random Forest Regressor. *Bioinformatics*. 2017;33(21):3477-3479.

9. Feng H et al. A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. *Nucleic Acids Res*. 2014;42(8):e69.

10. Wu H et al. Detection of differentially methylated regions from whole-genome bisulfite sequencing data without replicates. *Nucleic Acids Res*. 2015;43(21):e141.

11. Akalin A et al. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. *Genome Biol*. 2012;13:R87.

12. Hansen KD et al. BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions. *Genome Biol*. 2012;13:R83.

13. Hebestreit K et al. Detection of significantly differentially methylated regions in targeted bisulfite sequencing data. *Bioinformatics*. 2013;29(13):1647-1653.

14. Sun D et al. MOABS: model based analysis of bisulfite sequencing data. *Genome Biol*. 2014;15:R38.

15. Jühling F et al. metilene: fast and sensitive calling of differentially methylated regions from bisulfite sequencing data. *Genome Res*. 2016;26(2):256-262.

16. Yu X & Sun S. HMM-DM: identifying differentially methylated regions using a hidden Markov model. *Stat Appl Genet Mol Biol*. 2016;15(1):69-81.

17. Smyth GK. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Stat Appl Genet Mol Biol*. 2004;3:Article3.

18. Peters TJ et al. De novo identification of differentially methylated regions in the human genome. *Epigenetics Chromatin*. 2015;8:6.

19. Jaffe AE et al. Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies. *Int J Epidemiol*. 2012;41(1):200-209.

20. Phipson B et al. missMethyl: an R package for analyzing data from Illumina's HumanMethylation450 platform. *Bioinformatics*. 2016;32(2):286-288.

21. Park Y et al. Comprehensive Evaluation of Differential Methylation Analysis Methods for Bisulfite Sequencing Data. *Int J Environ Res Public Health*. 2021;18:7975.

22. Shafi A et al. A survey of the approaches for identifying differential methylation using bisulfite sequencing data. *Brief Bioinform*. 2018;19(5):737-753.

23. Gatev E et al. An evaluation of the genome-wide false positive rates of common methods for identifying differentially methylated regions using Illumina methylation arrays. *Epigenetics*. 2022.

24. Abramov S et al. Assessing the Differential Methylation Analysis Quality for Microarray and NGS Platforms. *Int J Mol Sci*. 2023;24:8591.

25. Maksimovic J et al. A cross-package Bioconductor workflow for analysing methylation array data. *F1000Research*. 2016;5:1281.

---

*Document prepared: December 2025*

*This report is intended as a technical reference for researchers selecting and implementing differential methylation analysis methods. Method selection should be based on data type, study design, and analytical goals.*
