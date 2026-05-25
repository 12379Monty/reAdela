# cfMeDIP-seq Protocol Reference

Shen et al. 2019 · Burgener et al. 2021 · Wilson et al. 2022  
4 Stages · 15 Steps · QC Metrics

---

## STAGE 1 — Library Preparation (Day 1 · Steps 1–3)

### Protocol Steps

| # | Step Name | Processing | Duration |
|---|-----------|------------|----------|
| 1 | End-Repair & A-Tailing | Input 1–10 ng plasma cfDNA. Enzymatic repair creates blunt ends with 5′ phosphates; T4 DNA polymerase and Klenow add 3′ adenosine (A) overhangs required for adapter ligation. | ~2–3 h |
| 2 | Adapter Ligation [UMI] | Ligate methylated Y-shaped sequencing adapters to the A-tailed library. In UMI implementations (Burgener et al. 2021), adapters carry an embedded 8–10 nt random barcode for downstream PCR deduplication. UMI option: barcode incorporated here; no extra wet-lab step required. | ~2–3 h |
| 3 | SPRI Bead Purification | Solid-Phase Reversible Immobilization (SPRI) cleanup removes excess adapters, enzymes, and small fragments. Bead-to-sample ratio controls size cutoff. | ~1 h |

### Sources of Variability

- **MODERATE** — Input DNA quantity and quality: low or degraded input amplifies stochasticity throughout the protocol
- **MODERATE** — Ligation efficiency: affected by adapter concentration, enzyme lot, and temperature; drives adapter dimer formation
- **LOW** — Enzymatic reaction temperature and timing (end-repair, ligation)
- **LOW** — SPRI bead-to-DNA ratio: affects recovery yield and adapter removal completeness
- **LOW** — Ethanol wash quality and evaporation timing during SPRI

### QC Metrics

| Metric | Description | Input Material / Method | Result Interpretation | QC Limit / Target |
|--------|-------------|-------------------------|-----------------------|-------------------|
| cfDNA input concentration | Starting cfDNA mass before Step 1; establishes baseline for all downstream yields | Pre-Step 1 plasma cfDNA · Qubit fluorometry | Low input (<1 ng) → stochastic library; inconsistent CpG coverage | ≥1 ng (target 1–10 ng) |
| Post-ligation yield | DNA mass recovered after adapter ligation and first SPRI cleanup | Post-Step 2–3 · Qubit | Low yield → insufficient library; may not produce adequate sequencing reads | ≥80% recovery of input |
| Adapter dimer content | Percentage of ~120–150 bp fragments arising from adapter self-ligation | Post-Step 3 · Bioanalyzer / TapeStation (optional) | High dimer content → wasted sequencing capacity; reduced effective read depth | <20% dimers; no dimer peak visible |

---

## STAGE 2 — Methylated DNA Immunoprecipitation (Day 2 · Steps 4–10)

⚠ Core of assay; highest variability

### Protocol Steps

| # | Step Name | Processing | Duration |
|---|-----------|------------|----------|
| 4 | Filler DNA + Spike-In Addition | Add filler DNA to bring total DNA to ~1 µg. Filler options: Lambda phage (unmethylated) or human genomic DNA (methylated). Add Arabidopsis thaliana spike-ins (methylated + unmethylated) or synthetic spike-ins (Wilson et al. 2022) for specificity calibration. Filler type is the single largest source of inter-lab batch effects. | ~30 min |
| 5 | DNA Denaturation | Heat to 95°C for 10 min to convert dsDNA to ssDNA accessible to anti-5mC antibody. Immediately chill on ice to prevent re-annealing. | ~15 min |
| 6 | Anti-5mC Immunoprecipitation | Incubate denatured DNA with 5 µg anti-5-methylcytosine (5mC) antibody overnight at 4°C with continuous rotation. Antibody selectively binds methylated ssDNA fragments. | ~16–18 h |
| 7 | Protein A/G Bead Capture | Add Protein A/G-conjugated magnetic beads to capture antibody–DNA complexes. Incubate 2 h at 4°C with rotation. | ~2 h |
| 8 | Washing | Wash beads 3–5× with IP wash buffer to remove non-specifically bound DNA and unbound antibody. Number and vigor of washes directly control background signal. | ~30 min |
| 9 | Elution | Elute captured DNA from beads via proteinase K digestion at 50°C for 2–3 h to degrade the antibody and release DNA. | ~3 h |
| 10 | Purification of IP DNA | SPRI bead cleanup of eluted DNA. Primary measurement checkpoint for Stage 2: IP recovery (Qubit) and specificity/fold-enrichment (qPCR) are measured here. | ~1 h |

### Sources of Variability

- **HIGH** — Filler DNA type (methylated vs unmethylated): drives systematic differences in IP efficiency; accounts for ~5% of variance (Wilson et al. 2022); major inter-lab batch effect
- **HIGH** — Antibody lot-to-lot variation: different lots show different binding efficiencies; affects all samples in a batch uniformly; difficult to correct without spike-in standards
- **MODERATE** — Wash stringency: number and vigor of washes controls background; inconsistency introduces run-to-run variability
- **MODERATE** — Antibody concentration and incubation conditions (temperature, rotation speed, mixing efficiency)
- **LOW** — Denaturation temperature accuracy and cooling rate
- **LOW** — Bead quality, non-specific binding, and magnetic separation efficiency
- **LOW** — Proteinase K activity and elution temperature/timing

### QC Metrics (all measured at Step 10, post-elution)

| Metric | Description | Input Material / Method | Result Interpretation | QC Limit / Target |
|--------|-------------|-------------------------|-----------------------|-------------------|
| IP recovery yield | Total DNA mass recovered after elution and purification; reflects overall IP efficiency | Post-Step 10 eluate · Qubit | Very low yield → IP failure; insufficient template for amplification | Sufficient for PCR (typically 0.1–5 ng) |
| Methylation specificity ▶ GO/NO-GO | % of captured DNA that is methylated, calculated from methylated vs unmethylated spike-in recovery | Post-Step 10 · qPCR with spike-in primers | <95% → poor IP enrichment; high unmethylated background; do not proceed | **>95%** |
| Fold-enrichment ratio ▶ GO/NO-GO | Ct(HIST1H2BA, methylated locus) relative to Ct(GAPDH, unmethylated locus); measures IP specificity for endogenous targets | Post-Step 10 · qPCR | <25-fold → IP failure or antibody lot problem; repeat IP or change antibody lot | **>25-fold** |
| Spike-in recovery balance | Ratio of reads/signal from methylated vs unmethylated spike-ins; batch effect diagnostic | Post-Step 10 · qPCR with spike-in primers | Imbalanced ratio → filler DNA or antibody batch effect; flags need for normalization | >19:1 (meth:unmeth) |

---

## STAGE 3 — Library Amplification (Day 3 · Steps 11–14)

### Protocol Steps

| # | Step Name | Processing | Duration |
|---|-----------|------------|----------|
| 11 | qPCR Specificity Check | Perform qPCR on a small aliquot of the IP eluate. Measure Ct for methylated locus (HIST1H2BA, positive control) and unmethylated locus (GAPDH, negative control). Calculate fold-enrichment and spike-in specificity. Formal Go/No-Go gate before committing to PCR amplification. | ~2 h |
| 12 | PCR Amplification [UMI] | Amplify IP library with indexed primers and a high-fidelity polymerase. Cycle number (typically 12–15) is determined empirically. Fewer cycles preserve library complexity; more cycles increase duplication rate and GC/length bias. UMI: barcodes incorporated at Step 2; post-sequencing deduplication (umi_tools dedup) removes PCR duplicates while retaining biological duplicates. | ~2 h |
| 13 | Final Purification & Size Selection | SPRI bead cleanup of PCR product. Perform size selection (200–600 bp) to remove residual adapter dimers (~120–150 bp) and very large fragments. | ~1 h |
| 14 | Library Quantification & QC | Quantify final library by Qubit or qPCR (nM). Assess fragment size distribution by Bioanalyzer or TapeStation. Check for adapter dimer peaks. Second formal Go/No-Go gate before sequencing. | ~1 h |

### Sources of Variability

- **MODERATE** — PCR cycle number: more cycles increase duplication rate, GC-content bias, and fragment length bias; systematic within a batch but variable between batches
- **MODERATE** — PCR GC bias: high-fidelity polymerases reduce but do not eliminate biased amplification of GC-rich methylated regions
- **LOW** — qPCR technical variability and primer efficiency: affects accuracy of Go/No-Go decision
- **LOW** — Size selection SPRI ratio: stringency directly affects fragment size profile and adapter dimer removal
- **LOW** — Quantification method accuracy: Qubit vs qPCR can differ; affects library loading for sequencing

### QC Metrics

| Metric | Description | Input Material / Method | Result Interpretation | QC Limit / Target |
|--------|-------------|-------------------------|-----------------------|-------------------|
| Ct value (qPCR) | PCR cycle threshold for HIST1H2BA and GAPDH amplification from IP aliquot; surrogate for IP recovery | IP eluate aliquot · Step 11 qPCR | High Ct (>30) → insufficient IP recovery; should not proceed without investigation | Ct <30 for methylated target |
| Fold-enrichment (qPCR gate) ▶ GO/NO-GO | HIST1H2BA vs GAPDH Ct ratio confirming IP specificity before amplification | IP eluate aliquot · Step 11 qPCR | <25-fold → do not proceed; re-do IP or change antibody lot | **>25-fold** |
| Specificity (qPCR gate) ▶ GO/NO-GO | Spike-in–based specificity confirming <5% unmethylated contamination | IP eluate aliquot · Step 11 qPCR | <95% → do not proceed; IP background too high | **>95%** |
| PCR cycles used | Number of amplification cycles applied; recorded as metadata for downstream correction | Step 12 — operator record | >15 cycles → elevated duplication rate and GC/length bias in sequencing data | 12–15 cycles |
| Fragment size peak | Modal fragment size of the amplified library on Bioanalyzer/TapeStation electropherogram | Post-Step 14 · Bioanalyzer or TapeStation | Peak outside 250–400 bp → size selection problem; adapter dimer carry-over | Peak 250–400 bp (range 100–700 bp) |
| Adapter dimer peak | Presence of ~120–150 bp peak on electropherogram | Post-Step 14 · Bioanalyzer or TapeStation | Present → wasted sequencing reads; repeat size-selection SPRI step | Absent (no peak at ~150 bp) |
| Library concentration ▶ GO/NO-GO | Final library molarity before sequencing submission | Post-Step 14 · Qubit or qPCR (nM) | <2 nM → insufficient for sequencing; re-amplify or repeat library prep | **≥2 nM** |
| UMI diversity [UMI] | Diversity of UMI barcodes in the library; low diversity indicates over-amplification or low-complexity library | Post-sequencing · umi_tools dedup output | Low diversity → most reads are PCR duplicates; library complexity insufficient | High diversity (most barcodes unique) |

---

## STAGE 4 — High-Throughput Sequencing (Day 4+ · Step 15)

### Protocol Steps

| # | Step Name | Processing | Duration |
|---|-----------|------------|----------|
| 15 | Paired-End Sequencing [UMI] | Pool indexed libraries by equimolar concentration. Sequence on Illumina platform (HiSeq, NovaSeq, or NextSeq) with paired-end reads (2×75 bp or 2×150 bp). Target depth 5–20 million read pairs per sample. Bioinformatic pipeline: FastQC → Trim Galore → (umi_tools extract if UMIs) → BWA-MEM alignment → SAMtools sort/filter → (umi_tools dedup if UMIs, else Picard MarkDuplicates) → MEDIPS or QSEA for CpG-level quantification. UMI: umi_tools extract moves barcode into read name before alignment; umi_tools dedup uses barcode + mapping coordinate to remove PCR duplicates while retaining true biological duplicates. | Day 4+ |

### Sources of Variability

- **HIGH** — Sequencing depth: lower depth reduces sensitivity for rare methylation events and low tumor-fraction samples; predictable effect that can be partially modeled
- **MODERATE** — Sequencing quality (Q-scores): base-calling errors introduce mapping artefacts; lane-to-lane and run-to-run variation on Illumina platforms
- **MODERATE** — Cluster density: over/under-clustering reduces quality scores and mapping rates
- **MODERATE** — Fragment length, GC content, CpG density: biophysical properties affect antibody capture (Stage 2) and are reflected in read distributions; can be modeled with normalization
- **LOW** — Library pooling imbalance: unequal molar concentrations lead to uneven depth per sample

### QC Metrics

| Metric | Description | Input Material / Method | Result Interpretation | QC Limit / Target |
|--------|-------------|-------------------------|-----------------------|-------------------|
| Raw read count | Total paired-end read pairs generated; sets the effective detection ceiling | FASTQ files · Sequencer output | Below target → insufficient sensitivity, especially for low tumor-fraction samples | 5–20 M paired-end reads |
| Q30 percentage | % of bases with Phred quality score ≥30 (≤0.1% error rate); reflects overall run quality | FASTQ · FASTQC | <85% → low-quality run; increased mapping errors and false-positive methylation calls | >85% |
| Adapter content | % of reads containing residual adapter sequences post-trimming | FASTQ · FASTQC | >5% → incomplete adapter trimming; biased alignment of short inserts | <5% |
| Mapping rate ▶ INCLUDE/EXCLUDE | % of reads that align to the reference genome (hg38) | BAM · BWA-MEM / SAMtools flagstat | <90% → contamination, low-quality library, or adapter issues; exclude or re-sequence | **>90%** |
| Duplication rate ▶ INCLUDE/EXCLUDE | % of reads flagged as PCR or optical duplicates | BAM · Picard MarkDuplicates or umi_tools dedup | >30% → excessive PCR amplification; data quality compromised; may warrant exclusion | **<30% (ideally <20%)** |
| Properly paired rate | % of read pairs with correct orientation and concordant alignment | BAM · SAMtools flagstat | <85% → adapter ligation problem, library structure issue, or contaminating sequences | >85% |
| Insert size distribution | Fragment length profile of sequenced inserts; should reflect nucleosomal cfDNA architecture | BAM · Picard CollectInsertSizeMetrics | Non-nucleosomal peak or very broad distribution → library degradation or size-selection failure | Peak ~166 bp (mononucleosomal); range 100–700 bp |
| CpG enrichment ▶ INCLUDE/EXCLUDE | Ratio of observed CpG frequency in reads vs expected genome-wide CpG frequency; measures methylation capture success | BAM · MEDIPS (CpG.enrichment) | <1.5 → poor methylation enrichment; antibody or IP failure reflected in sequencing data | **>1.5** |
| Saturation | Whether additional sequencing reads would capture new genomic regions; determines if depth is adequate | BAM · MEDIPS saturation analysis | Not at plateau → library under-sequenced; additional reads would improve coverage | Approaching saturation plateau |
| Spike-in read counts | Reads mapping to Arabidopsis thaliana or synthetic spike-in sequences; used for inter-sample normalization | BAM · BWA alignment to spike-in genome | Inconsistent counts across samples → batch effects from Stage 2 (filler type or antibody lot) | Consistent proportion across all samples in batch |
| Spike-in meth/unmeth ratio | Ratio of reads from methylated spike-ins to unmethylated spike-ins; directly mirrors IP specificity in sequencing data | BAM · Custom spike-in demultiplexing | <19:1 → poor IP specificity carried into sequencing data; flag sample for review | >19:1 (meth:unmeth) |
| UMI dedup rate [UMI] | % of reads removed by UMI-aware deduplication; distinguishes PCR duplicates from true biological duplicates | BAM · umi_tools dedup output log | Very high dedup rate (>70%) → low library complexity; may indicate IP failure or over-amplification | Consistent with PCR cycles used; <50% expected at 12–15 cycles |
