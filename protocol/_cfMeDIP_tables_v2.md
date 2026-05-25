<!--
# cfMeDIP-seq Protocol — Stage Tables
Shen et al. 2019 [@Shen:2018aa] · Burgener et al. 2021 [@Burgener:2021aa] · 
Wilson et al. 2022 [@Wilson:2022aa]
-->

#### STAGE 1 — Library Preparation (Day 1 · Steps 1–3)

##### Steps

| # | Name | Processing | Time |
|---|------|-----------|------|
| 1 | End-Repair & A-Tailing | Repair ends to blunt; add 3′ A-overhang; input 1–10 ng cfDNA | 2–3 h |
| 2 | Adapter Ligation | Ligate methylated adapters carrying unique molecular barcodes | 2–3 h |
| 3 | SPRI Purification | SPRI bead cleanup; remove excess adapters | 1 h |

##### Variability

| Level | Source |
|-------|--------|
| MODERATE | DNA input quantity and integrity |
| MODERATE | Ligation efficiency; adapter dimer formation |
| LOW | SPRI bead-to-DNA ratio; ethanol wash quality |

##### Metrics

| Metric | Description | Method | Interpretation | Limit |
|--------|-------------|--------|----------------|-------|
| cfDNA input | Starting mass | Qubit, pre-Step 1 | <1 ng → stochastic library | ≥1 ng |
| Post-ligation yield | Recovery after ligation + SPRI | Qubit, post-Step 3 | Low → insufficient library | ≥80% of input |
| Adapter dimer % | Self-ligated adapter fragments | Bioanalyzer, post-Step 3 | High → wasted reads | <20%; no dimer peak |

---

#### STAGE 2 — Immunoprecipitation (Day 2 · Steps 4–10)

##### Steps

| # | Name | Processing | Time |
|---|------|-----------|------|
| 4 | Filler DNA + Spike-Ins | Add filler to ~1 µg total; add methylated + unmethylated spike-in controls | 30 min |
| 5 | Denaturation | Heat to 95°C for 10 min → ssDNA; chill on ice | 15 min |
| 6 | Anti-5mC IP | Incubate with 5 µg anti-5mC antibody, overnight 4°C with rotation | 16–18 h |
| 7 | Bead Capture | Protein A/G beads capture antibody–DNA complex; 2 h, 4°C | 2 h |
| 8 | Washing | 3–5× IP buffer washes; remove non-specific DNA | 30 min |
| 9 | Elution | Proteinase K digestion at 50°C releases IP DNA | 2–3 h |
| 10 | IP Purification | SPRI cleanup of eluate; primary measurement checkpoint | 1 h |

##### Variability

| Level | Source |
|-------|--------|
| HIGH | Filler DNA type — dominant systematic batch effect (~5% variance; Wilson 2022) |
| HIGH | Antibody lot — batch-wide binding efficiency shift; requires spike-ins to detect |
| MODERATE | Wash stringency — controls background; inconsistency → run-to-run variation |
| MODERATE | Antibody concentration; incubation temperature and rotation speed |
| LOW | Denaturation accuracy; bead non-specific binding |

##### Metrics *(measured at Step 10, post-elution)*

| Metric | Description | Method | Interpretation | Limit |
|--------|-------------|--------|----------------|-------|
| IP recovery yield | DNA mass post-elution | Qubit | Very low → IP failure | ~0.1–5 ng |
| Methylation specificity ▶ GATE | % captured DNA that is methylated | qPCR, spike-in primers | <95% → poor enrichment; do not proceed | **>95%** |
| Fold-enrichment ▶ GATE | HIST1H2BA vs GAPDH Ct ratio | qPCR | <25 → IP or antibody failure | **>25-fold** |
| Spike-in balance | Meth vs unmeth spike-in signal | qPCR, spike-in primers | Imbalanced → batch effect present | >19:1 |

---

#### STAGE 3 — Amplification (Day 3 · Steps 11–14)

##### Steps

| # | Name | Processing | Time |
|---|------|-----------|------|
| 11 | qPCR QC | qPCR on IP aliquot: HIST1H2BA (positive) vs GAPDH (negative); Go/No-Go gate | 2 h |
| 12 | PCR Amplification | 12–15 cycles with indexed primers and high-fidelity polymerase | 2 h |
| 13 | Size Selection | SPRI cleanup; size-select 200–600 bp; remove adapter dimers | 1 h |
| 14 | Library QC | Quantify by Qubit/qPCR; assess size and dimers by Bioanalyzer/TapeStation | 1 h |

##### Variability

| Level | Source |
|-------|--------|
| MODERATE | PCR cycle number — more cycles → higher duplication rate and GC bias |
| MODERATE | GC bias — preferential amplification of GC-rich methylated regions |
| LOW | qPCR primer efficiency; SPRI size-selection ratio; quantification method |

##### Metrics

| Metric | Description | Method | Interpretation | Limit |
|--------|-------------|--------|----------------|-------|
| Ct value | Cycle threshold for HIST1H2BA | qPCR, Step 11 aliquot | Ct >30 → low IP yield | <30 |
| Fold-enrichment ▶ GATE | HIST1H2BA vs GAPDH | qPCR, Step 11 | <25 → do not amplify; repeat IP | **>25-fold** |
| Specificity ▶ GATE | Spike-in enrichment specificity | qPCR, Step 11 | <95% → do not amplify | **>95%** |
| PCR cycles | Cycles applied | Operator record | >15 cycles → elevated duplication and bias | 12–15 cycles |
| Fragment size peak | Modal fragment size | Bioanalyzer, Step 14 | Outside 250–400 bp → size-selection issue | 250–400 bp |
| Adapter dimer peak | ~150 bp electropherogram peak | Bioanalyzer, Step 14 | Present → repeat SPRI cleanup | Absent |
| Library concentration ▶ GATE | Final library molarity | Qubit/qPCR, Step 14 | <2 nM → re-amplify or repeat | **≥2 nM** |

---

#### STAGE 4 — Sequencing (Day 4+ · Step 15)

##### Steps

| # | Name | Processing | Time |
|---|------|-----------|------|
| 15 | Paired-End Sequencing | Pool equimolar libraries; Illumina PE 2×75 or 2×150 bp; 5–20 M read pairs per sample | Day 4+ |

##### Variability

| Level | Source |
|-------|--------|
| HIGH | Sequencing depth — directly limits sensitivity for low tumor-fraction samples |
| MODERATE | Base-call quality; cluster density; lane-to-lane run variation |
| MODERATE | Fragment length, GC, CpG density biases from Stage 2 reflected in read distribution |
| LOW | Library pooling imbalance |

##### Metrics

| Metric | Description | Method | Interpretation | Limit |
|--------|-------------|--------|----------------|-------|
| Raw read count | Total read pairs | Sequencer output | Below target → insufficient sensitivity | 5–20 M pairs |
| Q30 % | Bases with Phred ≥30 | FASTQC | <85% → low-quality run | >85% |
| Adapter content | Residual adapters post-trim | FASTQC | >5% → trim failure | <5% |
| Mapping rate ▶ GATE | % reads aligned to hg38 | BWA / SAMtools flagstat | <90% → contamination or library failure | **>90%** |
| Duplication rate ▶ GATE | % PCR / optical duplicates | Picard / umi_tools | >30% → over-amplification | **<30%** |
| Properly paired rate | % concordant read pairs | SAMtools flagstat | <85% → library structure issue | >85% |
| Insert size peak | Modal insert length | Picard | Off-nucleosomal → size-selection failure | ~166 bp |
| CpG enrichment ▶ GATE | CpG freq. in reads vs genome | MEDIPS | <1.5 → IP failure; exclude sample | **>1.5** |
| Saturation | Additional reads → new coverage? | MEDIPS | Not saturated → under-sequenced | Approaching plateau |
| Spike-in read counts | Reads mapping to spike-in sequences | BWA to spike-in genome | Inconsistent → Stage 2 batch effect | Consistent across batch |
| Spike-in meth/unmeth ratio | Meth vs unmeth spike-in reads | Custom analysis | <19:1 → poor IP specificity | >19:1 |

---
*▶ GATE = formal Go/No-Go decision point; bold limits are hard thresholds*
