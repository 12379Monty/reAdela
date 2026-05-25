<!--
## DNA Methylation Company Landscape {.tabset}
-->

* **DNA Methylation-Based Cancer Diagnostics: Company Landscape**

> Regulatory status as of early 2026.

**Abbreviations**

- **BSC** — bisulfite conversion (sodium bisulfite treatment converts unmethylated cytosines to uracil; leaves 5-methylcytosine intact)
- **IP** — immunoprecipitation (anti-5mC antibody pull-down enriches methylated fragments; no chemical conversion)
- **cfMeDIP-seq** — cell-free methylated DNA immunoprecipitation + high-throughput sequencing
- **MS-PCR** — methylation-specific PCR (bisulfite-converted DNA amplified with allele-specific primers)
- **NGS** — next-generation sequencing
- **WG** — whole-genome scope
- **LP** — large panel (>10K CpG sites; targeted enrichment)
- **SP** — small panel (<1K CpG sites or regions)
- **LDT** — laboratory-developed test (CLIA-certified lab, not FDA-cleared)
- **TOO** — tissue of origin
- **MRD** — minimal residual disease

### Home {-}

 <br/>


### Table 1 — Company Overview {-}

| Company / Product | Use Case | Sample | Cancers | Chemistry | Scope | Regulatory / Stage | Geography |
|---|---|---|---|---|---|---|---|
| **Adela** / Adela MCED + MRD | Multi-cancer detection + MRD | Blood | Multi (>12) | IP (cfMeDIP-seq) + WGS | WG (whole methylome) | Development / clinical validation | US, Canada |
| **GRAIL** / Galleri® | Multi-cancer screening | Blood | Multi (>50) | BSC + targeted NGS | LP (~1M CpGs) | LDT; FDA Breakthrough | US, EU (NHS) |
| **Exact Sciences** / CancerGuard™ | Multi-cancer detection + CRC | Blood | Multi | BSC + NGS + immunoassay | LP + protein panel | Clinical validation | US (Abbott acq. pending) |
| **Guardant Health** / Shield MCD | Multi-cancer detection | Blood | Multi (bladder, CRC, lung, pancreatic) | BSC + targeted NGS | LP | LDT; FDA Breakthrough | US |
| **BillionToOne** / Northstar Response® | Treatment response monitoring | Blood | Pan-cancer | BSC + targeted NGS + QCT counting | LP (>2,200 loci) | LDT commercial | US (Nasdaq: BLLN) |
| **Burning Rock** / OverC™ / ELSA-seq | Multi-cancer detection + MRD | Blood | Multi (9–22) | BSC + targeted NGS (ELSA-seq) | LP (161,984 CpGs) | CE mark; FDA Breakthrough | China, EU, US (Nasdaq: BNR) |
| **Singlera Genomics** / PanSeer™ | Multi-cancer detection | Blood | Multi (5 types) | BSC + targeted NGS (haplotype) | LP (12,000 sites) | Clinical validation | US, China |
| **Nucleix** / Bladder EpiCheck® | Bladder surveillance & detection | Urine | Bladder | BSC + MS-PCR | SP (15 regions) | FDA 510(k); CE mark | US, EU |
| **Epigenomics AG** / Epi proColon® / Epi proLung® / HCCBloodTest | CRC, lung, HCC screening | Blood | CRC, lung, HCC | BSC + MS-PCR | SP (1 gene each) | FDA approved; CE mark; restructuring | US, EU, China |
| **MDxHealth** / AssureMDx™ | Bladder cancer diagnosis | Urine | Bladder | BSC + MS-PCR + mutation PCR | SP (3 methyl + 3 mut. genes) | LDT commercial | US, EU |
| **Precision Epigenomics** / EpiSeek™ | Multi-cancer screening | Blood | Multi (>12) | BSC + MS-PCR | SP (proprietary panel) | LDT commercial | US |
| **EarlyDiagnostics** / CancerRadar™ | Multi-cancer detection | Blood | Multi | BSC + targeted NGS (presumed) | LP | Development | US |
| **Genomictree** / EarlyTect® | Bladder & CRC detection | Urine / Blood | Bladder, CRC | BSC + MS-PCR | SP (targeted loci) | Approved (Korea) | Korea, East Asia |

<br/>

### Table 2 — Assay Detail {-}

| Company / Product | Factor | Detail |
|---|---|---|
| **Adela** / Adela MCED + MRD | Enrichment method | IP: anti-5mC antibody pull-down (cfMeDIP-seq); no bisulfite conversion; preserves cfDNA integrity |
| | Detection | WGS of immunoprecipitated (methylated) cfDNA fragments; ML classifier for detection, TOO, and MRD quantification |
| | Key loci / regions | Whole methylome (genome-wide; no targeted panel required) |
| | Distinguishing features | Only commercial-stage platform using IP rather than BSC; avoids bisulfite DNA damage; profiles all methylated fragments unbiasedly; same wet-lab assay reused across MCED and MRD applications; CAMPERR study (n>5,000, prospective); licensed from Princess Margaret Cancer Centre (De Carvalho lab) |
| **GRAIL** / Galleri® | Enrichment method | BSC of plasma cfDNA; hybrid-capture enrichment of >100,000 differentially methylated regions |
| | Detection | Targeted NGS (Illumina); neural network ML classifier |
| | Key loci / regions | ~1 million CpG sites across >100K genomic regions |
| | Distinguishing features | WGBS used in discovery (CCGA substudy 1); refined to targeted panel for clinical scalability; TOO classification across >50 cancer types; trained on CCGA cohort (n=15,254) |
| **Exact Sciences** / CancerGuard™ | Enrichment method | BSC of cfDNA (methylation arm); amplicon sequencing (mutation arm); immunoassay (protein arm) |
| | Detection | Multi-analyte ML integrating methylation score + mutation score + protein score |
| | Key loci / regions | Methylation: undisclosed targeted panel; Mutations: cancer driver genes; Proteins: CA125, CEA, HGF, OPN + others (CancerSEEK lineage) |
| | Distinguishing features | Three independent biomarker classes combined; ASCEND-2 (n=11,000): 50.9% sensitivity, 98.5% specificity |
| **Guardant Health** / Shield MCD | Enrichment method | BSC of plasma cfDNA; targeted capture enrichment |
| | Detection | Targeted NGS; TOO ML classifier |
| | Key loci / regions | Undisclosed targeted methylation panel; loci informative for bladder, CRC, lung, and pancreatic cancer |
| | Distinguishing features | FDA Breakthrough Device Designation (Jun 2025); average-risk population ages 45+ |
| **BillionToOne** / Northstar Response® | Enrichment method | BSC of plasma cfDNA; targeted enrichment of cancer-specific hypermethylated loci |
| | Detection | Targeted NGS with proprietary QCT (Quantitative Counting Template) spike-in controls enabling single-molecule absolute counting |
| | Key loci / regions | >2,200 cancer-specific hypermethylated loci (core assay: >500 loci) |
| | Distinguishing features | Absolute tumour burden quantification (not relative VAF); tissue-naive — no biopsy required; LOD down to 0.01% tumour fraction; validated across 12 cancer types (*Sci Rep* 2025) |
| **Burning Rock** / OverC™ / ELSA-seq | Enrichment method | BSC of plasma cfDNA; targeted capture of 161,984 pre-selected CpG sites |
| | Detection | Deep targeted NGS (ELSA-seq); ML classifier for detection and TOO |
| | Key loci / regions | 161,984 CpG sites (custom panel from public + in-house methylome data); MRD product: CanCatch® |
| | Distinguishing features | THUNDER study (n=2,385): 98.9% specificity, 69.1% sensitivity, 91.7% top-2 TOO accuracy; included in China's 2024 Primary Liver Cancer clinical guidelines |
| **Singlera Genomics** / PanSeer™ | Enrichment method | BSC of plasma cfDNA; targeted sequencing of 600 pre-selected genomic regions |
| | Detection | Targeted NGS with methylation haplotype calling (co-methylation pattern analysis on single DNA strands) |
| | Key loci / regions | 12,000 CpG sites across 600 regions |
| | Distinguishing features | Analyses coordinated multi-CpG patterns on single strands rather than per-site averages; detection demonstrated up to 4 years pre-clinical diagnosis |
| **Nucleix** / Bladder EpiCheck® | Enrichment method | BSC of urinary cfDNA; no capture enrichment |
| | Detection | qPCR with methylation-specific primers; continuous numeric output score 0–100 |
| | Key loci / regions | 15 genomic regions (sequences proprietary) |
| | Distinguishing features | CE marked + FDA 510(k) cleared; sensitivity 91% high-grade disease, NPV 99%; objective score avoids subjective cytology interpretation |
| **Epigenomics AG** / Epi proColon® 2.0 | Enrichment method | BSC of plasma cfDNA; no capture enrichment |
| | Detection | Real-time MS-PCR targeting SEPT9 promoter methylation |
| | Key loci / regions | SEPT9 promoter |
| | Distinguishing features | First FDA-approved blood-based CRC screening test; company currently in financial restructuring |
| **Epigenomics AG** / Epi proLung® | Enrichment method | BSC of bronchial lavage cfDNA |
| | Detection | MS-PCR targeting SHOX2 methylation |
| | Key loci / regions | SHOX2 gene |
| | Distinguishing features | CE marked; workup aid for suspected lung malignancy |
| **Epigenomics AG** / HCCBloodTest | Enrichment method | BSC of plasma cfDNA |
| | Detection | NGS bisulfite sequencing targeting SEPT7 methylation |
| | Key loci / regions | SEPT7 gene |
| | Distinguishing features | Sensitivity 76.7%, specificity 64.1% in clinical study NCT03804593 |
| **MDxHealth** / AssureMDx™ | Enrichment method | BSC of exfoliated urinary cells; no capture enrichment |
| | Detection | qPCR: MS-PCR for methylation loci + allele-specific PCR for mutation loci |
| | Key loci / regions | Methylation: OTX1, ONECUT2, TWIST1; Mutations: FGFR3, TERT, HRAS |
| | Distinguishing features | Hybrid methylation + mutation assay; logistic regression model including patient age; sensitivity 97%, specificity 83% |
| **Precision Epigenomics** / EpiSeek™ | Enrichment method | BSC of blood cfDNA; no capture enrichment |
| | Detection | MS-PCR with proprietary software-derived marker panel |
| | Key loci / regions | Proprietary panel (undisclosed); validated across lung, breast, colon, prostate and others |
| | Distinguishing features | Self-pay / functional medicine distribution via TruDiagnostic (~40K clinicians); lower cost than NGS-based competitors |
| **EarlyDiagnostics** / CancerRadar™ | Enrichment method | BSC of blood cfDNA (presumed) |
| | Detection | Targeted NGS + ML (limited public detail) |
| | Key loci / regions | Undisclosed |
| | Distinguishing features | Early-stage; limited published validation data available |
| **Genomictree** / EarlyTect® Bladder | Enrichment method | BSC of urinary cfDNA |
| | Detection | qPCR targeting ONECUT2 methylation |
| | Key loci / regions | ONECUT2 |
| | Distinguishing features | Approved in Korea; designed for urine-based surveillance |
| **Genomictree** / EarlyTect® Colon | Enrichment method | BSC of blood or stool cfDNA |
| | Detection | qPCR targeting CRC-associated methylation loci |
| | Key loci / regions | Undisclosed CRC methylation panel |
| | Distinguishing features | Commercially available in East Asia |

---

### Key Observations {-}

#### Technique landscape

All commercially deployed methylation-based cancer diagnostics use **bisulfite conversion (BSC)** as the core chemistry — with one important exception: **Adela**, which uses IP-based cfMeDIP-seq and explicitly avoids bisulfite conversion to prevent cfDNA degradation. Among the BSC-based players, the split is:

- **BSC + targeted NGS** (GRAIL, Burning Rock, Singlera, Guardant, Exact Sciences, BillionToOne): large panels captured by hybrid enrichment and sequenced; ML classifiers interpret the methylation signal. Most information-rich but highest cost. Best suited for multi-cancer detection and TOO localisation.
- **BSC + MS-PCR / qPCR** (Epigenomics, Nucleix, MDxHealth, Genomictree, Precision Epigenomics): bisulfite-treated DNA amplified with methylation-sensitive primers at a small number of pre-validated loci. Lower cost, faster turnaround, further along in regulatory clearance — but limited to single-cancer or narrow-scope applications.
- **BSC + multi-analyte** (Exact Sciences CancerGuard, MDxHealth AssureMDx): methylation combined with somatic mutation and/or protein biomarkers; ML integration boosts sensitivity at low tumour fractions.

#### Adela's IP-based distinction

Adela is the only company at commercial/clinical-validation stage using immunoprecipitation rather than bisulfite chemistry. The cfMeDIP-seq platform enriches methylated cfDNA fragments via an anti-5mC antibody pull-down before WGS, preserving far more input material than bisulfite conversion (which destroys ~96% of DNA). A key operational advantage is that the same wet-lab assay is reused unchanged across MCED and MRD applications — only the downstream ML classifier changes — enabling more efficient product development across the cancer care continuum.

#### Alternative chemistries — current status

| Chemistry | Principle | Why not commercial yet |
|---|---|---|
| **EM-seq** | Enzymatic conversion (TET2 + APOBEC3A) preserves DNA integrity better than BSC | NEB research kit available; no clinical validation datasets; no Dx company has committed |
| **Nanopore direct sequencing** | Native 5mC detection from ionic current; no conversion; long reads enable methylation haplotypes | cfDNA input sensitivity not yet validated for clinical use |
| **Illumina EPIC / 850K array** | BSC + bead array; ~850K CpGs | Requires micrograms of input DNA; not suitable for cfDNA liquid biopsy |
| **Restriction enzyme-based (MSRE, IMPRESS)** | Methylation-sensitive enzymes cut unmethylated sites; bisulfite-free | Academic proof-of-concept only; not yet clinical |

#### Genomic scope vs. commercial maturity

The smaller the panel, the more commercially mature: single-gene MS-PCR tests (SEPT9, SHOX2) achieved FDA approval first; large-panel NGS tests are still operating as LDTs or navigating the PMA pathway. The trade-off is sensitivity and cancer-type breadth — a large-panel or genome-wide approach is required to detect low-abundance early-stage methylation signals across multiple tumour types simultaneously.

---

*Table compiled from public sources. Last updated: April 2026.*
