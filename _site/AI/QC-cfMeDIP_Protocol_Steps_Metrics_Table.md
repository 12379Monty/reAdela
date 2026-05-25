# cfMeDIP-seq Protocol: Steps, Stages, and Derived Metrics

## Protocol Steps and Quality Metrics Table

| Step | Stage | Step Name | Primary Metrics Derived | Measurement Method | Typical Threshold/Target |
|------|-------|-----------|------------------------|--------------------|-----------------------|
| 1 | 1: Library Prep | End-repair & A-tailing | DNA concentration post-reaction | Qubit fluorometry | ≥80% recovery of input |
| 2 | 1: Library Prep | Adapter ligation | Ligation efficiency; adapter-dimer content | qPCR or Bioanalyzer | >50% ligation efficiency; <20% dimers |
| 3 | 1: Library Prep | SPRI bead purification | DNA recovery yield; residual adapter content | Qubit; Bioanalyzer | ≥70% recovery; no adapter peak |
| 4 | 2: Immunoprecipitation | Add filler DNA & spike-ins | Total DNA mass; spike-in concentration | Qubit | ~1 μg total DNA; defined spike-in amount |
| 5 | 2: Immunoprecipitation | Denaturation (95°C) | *(No direct metric - process step)* | — | — |
| 6 | 2: Immunoprecipitation | Anti-5mC IP (overnight) | *(Metrics derived post-elution)* | — | — |
| 7 | 2: Immunoprecipitation | Bead capture & wash | *(Metrics derived post-elution)* | — | — |
| 8 | 2: Immunoprecipitation | Elution & purification | IP DNA recovery; methylation specificity; fold-enrichment ratio | Qubit; qPCR | Specificity >95%; fold-enrichment >25 |
| 9 | 3: Amplification | qPCR QC | Methylation specificity; fold-enrichment (HIST1H2BA/GAPDH); Ct values | qPCR | Specificity >95%; fold-enrichment >25; Ct <30 |
| 10 | 3: Amplification | PCR amplification | Amplified library yield; cycle number used | Qubit | Sufficient yield for sequencing; typically 12-15 cycles |
| 11 | 3: Amplification | Size selection & QC | Fragment size distribution; library concentration; library molarity | Bioanalyzer/TapeStation; Qubit/qPCR | Peak 250-400 bp; no adapter dimers; ≥2 nM |
| 12 | 4: Sequencing | Paired-end sequencing | See detailed sequencing metrics below | Sequencer output; FASTQC; MEDIPS | See below |

---

## Detailed Metrics by Stage

### Stage 1: Library Preparation Metrics

| Metric | Description | Source Step | Concern if Abnormal |
|--------|-------------|-------------|---------------------|
| Input DNA concentration | Starting cfDNA amount | Pre-Step 1 | Low input → noisy data |
| Post-ligation yield | DNA recovered after adapter ligation | Step 2-3 | Low yield → insufficient library |
| Adapter dimer content | Small fragments from adapter self-ligation | Step 3 | High dimers → wasted sequencing reads |

### Stage 2: Immunoprecipitation Metrics

| Metric | Description | Source Step | Concern if Abnormal |
|--------|-------------|-------------|---------------------|
| Methylation specificity | % of captured DNA that is methylated (from spike-ins) | Step 8-9 | <95% → poor enrichment, high background |
| Fold-enrichment ratio | Ratio of methylated (HIST1H2BA) to unmethylated (GAPDH) signal | Step 8-9 | <25 → IP failure or antibody issue |
| IP recovery | Amount of DNA recovered post-IP | Step 8 | Very low → IP failure |
| Spike-in recovery | Recovery of methylated vs unmethylated spike-ins | Step 8-9 | Imbalanced → batch effect indicator |

### Stage 3: Amplification Metrics

| Metric | Description | Source Step | Concern if Abnormal |
|--------|-------------|-------------|---------------------|
| Ct value (qPCR) | Cycle threshold for amplification | Step 9 | High Ct (>30) → low IP recovery |
| PCR cycles used | Number of amplification cycles | Step 10 | >15 cycles → high duplication rate |
| Library concentration | Final library yield (ng/μL or nM) | Step 11 | Low → may need re-amplification |
| Fragment size peak | Mode of fragment size distribution | Step 11 | Outside 250-400 bp → size selection issue |
| Adapter dimer peak | Presence of ~120-150 bp peak | Step 11 | Present → wasted reads, repeat cleanup |

### Stage 4: Sequencing Metrics

| Metric | Description | Derived From | Typical Target |
|--------|-------------|--------------|----------------|
| **Raw read count** | Total read pairs generated | FASTQ files | 5-20 million pairs |
| **Q30 percentage** | % bases with quality score ≥30 | FASTQC | >85% |
| **Adapter content** | Residual adapter sequences | FASTQC | <5% |
| **Mapping rate** | % reads aligned to reference | BWA/SAMtools | >90% |
| **Duplication rate** | % duplicate reads (PCR or optical) | SAMtools/Picard | <30% (ideally <20%) |
| **Properly paired rate** | % read pairs correctly oriented | SAMtools | >85% |
| **Insert size distribution** | Fragment length profile | Picard | Peak ~166 bp (mononucleosomal) |
| **CpG enrichment** | Ratio of CpG frequency in reads vs genome | MEDIPS | >1.5 |
| **Saturation** | Whether additional reads would add information | MEDIPS | Approaching plateau |
| **Spike-in read counts** | Reads mapping to spike-in sequences | BWA alignment | Consistent across samples |
| **Spike-in methylated/unmethylated ratio** | Ratio of methylated to unmethylated spike-in reads | Custom analysis | >19:1 (reflects specificity) |

---

## Summary: Key Go/No-Go Decision Points

| Decision Point | Stage | Critical Metrics | Action if Failed |
|----------------|-------|------------------|------------------|
| Post-IP QC | 2→3 | Specificity >95%; Fold-enrichment >25 | Do not proceed; repeat IP |
| Pre-sequencing QC | 3→4 | Library concentration ≥2 nM; correct size profile | Re-amplify or repeat library prep |
| Post-sequencing QC | 4 | Mapping >90%; Duplication <30%; CpG enrichment >1.5 | Exclude sample or re-sequence |

---

## Metrics Workflow Diagram

```
Step 1-3 (Library Prep)
    │
    ├──→ DNA yield (Qubit)
    └──→ Size profile (Bioanalyzer) [optional at this stage]
    
Step 4-8 (IP)
    │
    └──→ Post-IP metrics measured at Step 8-9:
         ├── Specificity (spike-in qPCR)
         ├── Fold-enrichment (HIST1H2BA/GAPDH qPCR)
         └── Recovery yield (Qubit)
         
         ⚠️ GO/NO-GO DECISION
    
Step 9-11 (Amplification)
    │
    ├──→ Ct values (qPCR) [Step 9]
    ├──→ Library yield (Qubit) [Step 11]
    └──→ Size distribution (Bioanalyzer/TapeStation) [Step 11]
    
         ⚠️ GO/NO-GO DECISION
    
Step 12 (Sequencing)
    │
    └──→ Sequencing QC metrics:
         ├── Read quality (FASTQC)
         ├── Alignment stats (SAMtools)
         ├── Duplication (Picard)
         ├── Insert size (Picard)
         ├── CpG enrichment (MEDIPS)
         ├── Saturation (MEDIPS)
         └── Spike-in analysis
         
         ⚠️ INCLUDE/EXCLUDE DECISION
```

---

*This table supplements the cfMeDIP-seq Protocol, Pipeline, and Classifiers document.*
