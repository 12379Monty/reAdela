* Technical Review: Zeng et al. (2026)
   -  *A Pan-Cancer Compendium of 1,294 Plasma Cell-Free DNA Methylomes and Fragmentomes Enabling Multicancer Detection*
      -  Nature Cancer 7, 384–398

<!--
**Perspective:** Experienced statistician focused on data processing, normalization methodology, quantification approaches, and the development/validation structure of the analytical pipeline.
-->

---

### 1. Study Context and Motivation

This paper represents the field's first large-scale harmonization effort for cfMeDIP-seq data, aggregating 1,074 plasma cfMeDIP-seq profiles from nine distinct studies within The Cancer Genetics and Epigenetics (TCGE) project, spanning 11 cancer types, Li-Fraumeni syndrome (LFS) carriers, and healthy controls. The stated motivation is that broader application of cfMeDIP-seq has been limited by small cohorts and inconsistent data processing — a problem this paper directly addresses by developing a uniform computational workflow. Validation was performed in an independent set of 220 samples.

The paper is best understood as a **data resource and methods harmonization paper** rather than a biomarker validation study. Its primary technical contribution is the systematic evaluation and selection of a normalization strategy across heterogeneous cfMeDIP-seq datasets, and the characterization of pan-cancer methylome and fragmentome signatures from the resulting harmonized data.

---

### 2. Data Collection and Quality Control

#### 2.1 Dataset Composition

The study collated 1,074 cfMeDIP-seq profiles across nine studies comprising cancer samples from 11 cancer types, LFS carriers, and healthy controls.

The nine constituent studies introduce substantial technical heterogeneity:
- Sequencing modes: one study (TCGE-CFMe-MCA) used single-end (SE) sequencing; all others used paired-end (PE) sequencing
- Sequencing platforms and read lengths varied across studies
- Library preparation protocols, antibody lots, and filler DNA compositions differed across contributing labs

This heterogeneity is the central technical challenge the normalization strategy must address.

#### 2.2 Quality Control

Libraries were sequenced with a median of ~61 million raw sequencing reads and about 28 million unique reads after QC trimming and duplication removal. All samples were processed through MEDIPIPE v.1.1.0, which computed 21 comprehensive QC metrics including raw and preprocessed read depths, saturation, coverage, specificity, and enrichment scores.

More than 90% of the samples reached enrichment scores of GoGe greater than 1.7, relH greater than 3.0, and a saturation score (maxEstCor) exceeding 0.9 — thresholds previously suggested for high-quality data by the inventors of the cfMeDIP-seq method.

After QC, 974 high-quality samples (SE: 378; PE: 596) were retained for downstream analysis.

**Key observation:** A negative correlation was observed between enrichment scores (GoGe and relH) and fragment size in PE profiles, indicating that longer fragments may result in lower cfMeDIP signal-to-noise ratios. This is an important finding for normalization: fragment length affects not only IP efficiency (as Wilson et al. 2022 demonstrated) but also the overall enrichment quality metric used for QC gating.

---

### 3. Normalization and Quantification Methods

This is the section of greatest statistical interest. The paper systematically evaluated six normalization strategies before selecting one for all downstream analyses.

#### 3.1 The Six Candidate Normalization Strategies

The paper compared the following approaches, assessed primarily by their ability to reduce batch effects within the TCGE-CFMe-AML study (which included three technical replicates per participant across labs) and across four healthy control cohorts:

| Strategy | Description |
|----------|-------------|
| **Raw read count** | No normalization beyond alignment |
| **RPKM/FPKM** | Reads per kilobase per million; depth-normalized |
| **MEDEStrand estimated RMS** | Stranded sigmoid model for absolute methylation (MeDEStrand v.0.0.0.9000) |
| **QSEA estimated beta** | Sigmoidal model with empirical knowledge (QSEA v.1.20.0) |
| **DESeq2 normalized count** | Variance-stabilizing normalization without prior batch correction |
| **ComBat-seq + DESeq2** | Batch correction of raw counts followed by DESeq2 normalization |

**Selection criterion:** PCA-based assessment of batch variable association with top principal components, evaluated within-study (AML technical replicates) and across studies (healthy controls).

#### 3.2 Selected Strategy: ComBat-seq + DESeq2

Among six tested DNA methylation quantification and normalization strategies, the ComBat-seq + DESeq2 method most effectively reduced batch effects, especially when SE and PE samples were processed separately.

The workflow is:
1. **ComBat-seq** (sva v.3.52.0) applied to raw read counts, correcting labeled library preparation or sequencing batches while preserving sample subtype differences
2. **DESeq2** (v.1.44.0) normalization of the ComBat-seq corrected counts
3. Applied **separately** to SE and PE samples, since SE and PE samples showed intrinsic differences beyond sequencing batch

**Critical statistical observations:**

**1. This is a count-based normalization framework borrowed from RNA-seq.** ComBat-seq is designed for negative-binomial count data (RNA-seq read counts); DESeq2 uses median-of-ratios normalization, also designed for RNA-seq. Applying these tools to cfMeDIP-seq read counts treats enrichment-based methylation data as if it were expression count data. This is a pragmatic approximation — the NegBinom count model is plausible for binned cfMeDIP-seq read counts — but the biological meaning of DESeq2's size factors applied to immunoprecipitation-enriched data differs from the RNA-seq context for which they were designed. No mechanistic justification for this choice is provided beyond empirical batch reduction performance.

**2. The Stutheit-Zhao (2024) NegBinom EM mixture model for methylation probability is not used.** The output feature used for all downstream analyses is ComBat-seq + DESeq2 normalized read count — not methylation probability, not CSM, not picomolar molar amount. The paper thus works with normalized but relative count data, not absolute methylation estimates.

**3. The Wilson et al. (2022) spike-in normalization is not used**, despite Wilson being a co-author and the spike-in method being cited (reference 16). Spike-in normalization requires that spike-in DNA was added to all samples at library preparation, which was not the case for the historical datasets collected from nine independent studies — the spike-in approach is only applicable prospectively.

**4. MEDEStrand and QSEA are present in MEDIPIPE** but were not selected. Notably, MEDEStrand is used for a specific secondary purpose: PBL depletion filtering, where the median MEDEStrand-estimated absolute methylation value across 20 healthy PBL samples was used to identify PBL-depleted bins (those with median rms < 0.1). This is a methodologically interesting split — MEDEStrand estimates are used for the PBL filter but not for the primary quantification, presumably because MEDEStrand performed less well in the batch correction comparison.

**5. QSEA introduced spurious filler-associated variance.** The paper's Extended Data Fig. 2 (referenced in the Methods) shows the six-way comparison. This is consistent with the Wilson et al. 2022 finding that QSEA normalization made filler DNA composition significantly associated with PC1, a result that appeared to amplify rather than correct this technical variable.

#### 3.3 Bin Filtering

After normalization, bins were filtered to retain those:
- On autosomes only (sex chromosomes excluded)
- Not overlapping ENCODE blacklist regions
- Containing at least one CpG site
- For DMR identification: containing more than five CpG sites (n = 1,279,093 bins)

The final analysis space was **7,445,098 bins** containing at least one CpG site for the methylation quantification PCA, reduced to 1,279,093 bins with >5 CpGs for DMR identification. The top 10,000 variable bins (selected by interquartile range) were used for PCA visualization.

#### 3.4 What Is Not Corrected For

The ComBat-seq + DESeq2 approach explicitly corrects for labeled batch variables (library preparation batch, sequencing platform differences). It does **not** correct for:

- Fragment length biases in IP efficiency (addressed in Wilson et al. 2022 but not applicable here)
- GC content biases in enrichment
- CpG fraction effects on enrichment
- Sequencing depth differences beyond DESeq2's size factor

The approach relies on ComBat-seq to absorb these unlabeled sources of variation into its batch correction model, which is a strong assumption when batches are confounded with cancer type or study design.

---

### 4. DMR Identification: From Normalized Counts to Biological Features

#### 4.1 Statistical Framework

DMRs were identified using DESeq2 on ComBat-seq corrected data, requiring an absolute log₂ fold change >1 and a false discovery rate (FDR) <0.05.

The use of DESeq2 for both normalization and differential testing means that the same model assumptions underlie both steps. DESeq2's differential expression framework fits a NegBinom GLM per bin, shrinks dispersion estimates across bins using empirical Bayes, and applies the Wald test. This is a reasonable but not universally accepted approach for enrichment-based sequencing data.

#### 4.2 Filtering Pipeline for DMR Identification

DMR identification incorporated a multi-step filtering pipeline to remove non-cancer signals:

**Step 1 — Age and sex confounders:** DMRs identified without adjustment showed minimal overlap with age- (an average of 0.48% of hyper-DMRs and 0.40% of hypo-DMRs) and sex-associated (an average of 0.13% of hyper-DMRs and 0.17% of hypo-DMRs) signatures. The primary analysis did not include age and sex as DESeq2 covariates (due to missing data for >240 samples), but identified DMRs were filtered post-hoc against age- and sex-associated signatures.

**Step 2 — PBL depletion:** PBL-depleted bins (2,128,687 out of 7,445,098) were defined based on the criterion that the median MEDEStrand-estimated absolute methylation value across healthy PBLs was less than 0.1. Restricting to PBL-depleted bins retained on average 22.9% of hyper-DMRs and 19.4% of hypo-DMRs per cancer type — a substantial reduction that highlights the magnitude of immune cell background signal in cfMeDIP-seq.

**Step 3 — Blood and endothelial cell exclusion:** DMRs overlapping blood- or endothelial cell-associated regions were excluded: 0.52% of hyper-DMRs and 0.88% of hypo-DMRs.

**Step 4 — SE/PE intersection for the pan-cancer signature:** The 14,202 pan-cancer hyper-DMRs represent the intersection of SE and PE study DMRs. Of these, 99.8% appeared in multiple cancer types, and 66.4% in at least four.

#### 4.3 The Pan-Cancer Signature vs. Stutheit-Zhao (2024)

This pan-cancer signature of 14,202 common hyper-DMRs derived from cfMeDIP-seq data itself contrasts directly with the Stutheit-Zhao approach, which derived its 200-DMC signature from TCGA 450K array data. Zeng et al.'s signature is larger by two orders of magnitude (14,202 vs. 200 regions), derived from the same data type (cfMeDIP-seq), and identified through a formal differential testing framework rather than a TCGA array-based selection. The tradeoff is circularity: the signature is derived from the same data used for classification, unlike the fully external TCGA-based approach of Stutheit-Zhao.

---

### 5. Fragmentomic Features and Quantification

#### 5.1 Features Computed

For 473 PE baseline samples, four fragmentomic features were computed:

| Feature | Computation | Reference Method |
|---------|-------------|-----------------|
| **Fragment insert size** | Picard CollectInsertSizeMetrics; per-bin z-scores vs. healthy; NMF for latent signatures | Vessies et al. 2022 NMF approach |
| **Genome-wide fragment ratios** | Short (90–150 bp) to long (151–220 bp) within 5 Mb windows; GC and depth corrected | DELFI (Cristiano et al. 2019) |
| **Nucleosome footprinting** | Distance of 167 bp fragments from ~13M known nucleosome positions; z-scores vs. healthy | Vanderstichele et al. 2022 |
| **5′ end motifs** | Frequency of all 256 4-nucleotide sequences at fragment 5′ ends; z-scores vs. healthy | Jiang et al. 2020 |

All four features use **z-scoring against healthy controls** as their primary normalization/quantification approach, rather than absolute measurement. Per-sample genome-wide deviation is summarized as Σz (sum of per-bin or per-motif z-scores). This is a relative measure — the implicit assumption is that healthy control samples are stable and representative references.

#### 5.2 Fragment Score (FS) from NMF

The NMF approach factorizes insert size distributions into two latent components (cancer-enriched vs. healthy-enriched), computing a weighted fragment score measuring similarity to the cancer-associated profile. This is conceptually similar to the FLS in Stutheit-Zhao 2024, but computed genome-wide rather than restricted to signature windows. The NMF was trained on a 70% split of the data, with no formal held-out test for the NMF itself — the signature shapes are used descriptively.

#### 5.3 Correlation Among Fragmentomic Features

Features were partially correlated: FS tracked the short-fragment proportion (Spearman ρ = 0.89), and nucleosome peak-distance z-scores correlated with fragment ratios (ρ = 0.73). This correlation structure means that "integrating" all four features does not add four independent pieces of information. The PCA-based dimensionality reduction applied before classifier training partially addresses this, but the features are not independent inputs.

---

### 6. Development and Validation Structure

#### 6.1 Training/Test Framework for Classification

Classifiers were trained on principal components of cfDNA methylation and fragmentomic features, benchmarking seven algorithms by nested cross-validation, and tested on predefined independent datasets.

The classification framework uses:
- **Dimensionality reduction:** PCA applied within training folds, transformation applied to test folds to prevent data leakage
- **Algorithm benchmarking:** LASSO, SVM, random forest, gradient boosting, XGBoost, LDA, K-nearest neighbors
- **Cross-validation:** Nested k-fold, with fold counts scaled to dataset size (LOOCV for n<10; 5×5 for n=10–100; 5×10 for n=101–250; 10×10 for n>250)
- **Class imbalance:** Addressed by random down-sampling within resamples

#### 6.2 The Training/Validation Split

| Component | Training Data | Validation Data |
|-----------|--------------|-----------------|
| Pan-cancer methylation signature (14,202 DMRs) | All 9 TCGE studies (PE and SE separately) | ✗ Used in training |
| ComBat-seq + DESeq2 normalization | Fitted within each study | Applied to validation cohort |
| PCA transformation | Fitted on training fold | Applied to validation fold |
| Classifiers | Primary PE dataset (nested CV) | TCGE-CFMe-INSPIRE + TCGE-CFMe-HCC (n=95 baseline after QC) |

**Critical validation observation:** The combination of fragment ratios and methylation achieved the highest mAUC of 0.954 in the validation dataset, close to the primary PE results of 0.961. However, the most striking finding is the performance divergence for specific feature types in the validation set. End motifs achieved mAUC = 0.930 in the primary PE dataset but only 0.188 in the validation dataset — a catastrophic generalization failure. This is attributed to a different laboratory protocol for 17 healthy controls in the HCC study affecting 5′ end-motif distributions. This result powerfully illustrates how fragmentomic features, unlike methylation, are highly sensitive to pre-analytical and protocol variables.

**The INSPIRE cohort appears in the validation set** — which is the same cohort used as the primary data in Stutheit-Zhao 2024. Its TCGE-CFMe-INSPIRE data (EGAD00001011312) is used here as an independent validation cohort for the Zeng 2026 classifiers. The two papers thus use the INSPIRE data in opposite roles: as primary training/analysis data in Stutheit-Zhao and as a validation cohort in Zeng.

#### 6.3 Sample Size and Power

Power analyses were conducted, and the paper is unusually transparent in reporting them. For pan-cancer DMR identification in PE and SE samples, a moderate power of 73.7% was achieved with the sample sizes of 378 (SE) and 420 (PE). This means approximately 26% of true DMRs may have been missed at these sample sizes — a sobering figure for a 14,202-DMR signature. Integrative classifiers achieved 80% power with as few as 15–20 samples per group, reflecting the high dimensionality reduction achieved by PCA.

#### 6.4 The SE/PE Separation Problem

A major structural limitation that the authors handle explicitly but incompletely is the fundamental difference between SE and PE data. SE and PE samples showed intrinsic differences beyond sequencing, as confirmed by mimicked and original SE samples comparison, which also influenced differentially methylated regions identification. The authors applied ComBat-seq to correct SE/PE labels in addition to study labels, but the resulting pan-cancer signature (14,202 regions) relies on the intersection of SE and PE DMRs as a conservative filter. This conservatism is appropriate but means the pan-cancer signature may exclude DMRs detectable only in PE data (e.g., those dependent on fragment-level information).

---

### 7. Key Results in Quantitative Terms

| Result | Metric | Caveat |
|--------|--------|--------|
| Pan-cancer methylation classifier (PE, primary) | mAUC = 0.986 | Within-dataset nested CV |
| Combined methylation + fragment ratios (PE, primary) | mAUC = 0.961 | Within-dataset nested CV |
| Combined methylation + fragment ratios (PE, validation) | mAUC = 0.954 | Independent INSPIRE + HCC cohort |
| Methylation alone (PE, validation) | mAUC = 0.674 | Markedly lower than primary |
| End motifs (PE, validation) | mAUC = 0.188 | Protocol-driven failure to generalize |
| Sensitivity at 99% specificity (combined, validation) | 76.6 ± 4.0% | vs. 94.9% in primary dataset |

The drop from methylation mAUC 0.986 (primary) to 0.674 (validation) warrants careful attention. The validation cohort includes cancer types (melanoma, ovarian, breast, mixed) and healthy controls (HCC study, different protocol) not all represented in the training data — this is a legitimate test of generalization across cancer types and protocols.

---

### 8. Data and Code Availability

| Resource | Location |
|----------|----------|
| Source data (DMRs, PCA/UMAP plots, BED files) | Zenodo: https://doi.org/10.5281/zenodo.15191455 |
| All analysis code | GitHub: https://github.com/HansenHeLab/cfMeDIP-seq_Data_Resource_Codes |
| MEDIPIPE pipeline | https://github.com/yzeng-lol/MEDIPIPE |
| TCGE-CFMe-AML raw data | EGA: EGAS00001005069 |
| TCGE-CFMe-PRAD raw data | EGA: EGAS00001005522 |
| TCGE-CFMe-UM raw data | EGA: EGAD00001008998 |
| TCGE-CFMe-HBC + LFS raw data | EGA: EGAS00001006539 |
| TCGE-CFMe-INSPIRE raw data | EGA: EGAD00001011312 |
| TCGE-CFMe-HCC raw data | EGA: EGAD50000000652 |
| TCGE-CFMe-MCA, BCA, HNSC, SCLC | Upon request from corresponding authors |

The GitHub code repository and Zenodo deposit make the DMR BED files and source data for figures publicly accessible — an important resource for the field.

---

### 9. Summary of Key Statistical Considerations

| Issue | Severity | Detail |
|-------|----------|--------|
| RNA-seq normalization tools applied to enrichment-based methylation data | Moderate | ComBat-seq and DESeq2 were designed for RNA-seq count data; mechanistic justification for their use with cfMeDIP-seq is empirical, not principled |
| Pan-cancer DMR signature derived from same data used for classification | High | Circular use: 14,202 DMRs identified from TCGE studies are then used as features for classification in those same studies; only the within-fold PCA/classifier pipeline is held out |
| SE/PE intrinsic difference not fully correctable | Moderate | Two rounds of ComBat-seq partially address SE/PE differences but residual confounding likely remains; SE and PE classifiers are evaluated separately |
| End motif catastrophic generalization failure (mAUC 0.93 → 0.19) | High | Demonstrates extreme sensitivity of fragmentomic features to pre-analytical protocol differences; limits clinical applicability |
| Methylation generalization drop (mAUC 0.986 → 0.674) | Moderate | Substantial but expected given validation includes new cancer types and different protocols; combined classifier recovers to 0.954 |
| Power of 73.7% for DMR identification | Moderate | ~26% of true DMRs may be missed; the 14,202-DMR signature is already conservative |
| Missing data for age/sex covariates (>240 samples) | Low–Moderate | Age and sex not included as DESeq2 covariates in primary analysis; filtered post-hoc using external signatures, which is a reasonable but imperfect substitute |
| Validation cohort overlap with prior publication | Low | INSPIRE cohort used as validation here was the primary cohort in Stutheit-Zhao 2024; not truly independent of the research group |

---

### 10. Relationship to Stutheit-Zhao (2024) and Wilson (2022)

This paper explicitly cites both Stutheit-Zhao 2024 (reference 14) and Wilson 2022 (reference 16), but adopts neither of their normalization approaches. The normalization lineage comparison across the three papers is:

| Paper | Methylation Quantification | Normalization Approach | Batch Correction |
|-------|---------------------------|----------------------|-----------------|
| Stutheit-Zhao 2024 | NegBinom EM mixture → methylation probability | Arithmetic sum (CSM); no explicit depth normalization | None (single-cohort) |
| Wilson 2022 | Gaussian GLM on spike-in read counts → picomoles | Absolute molar amount; corrects length, GC, CpG fraction | Spike-in-based |
| **Zeng 2026** | **ComBat-seq + DESeq2 normalized read counts** | **Relative count; corrects labeled batches** | **ComBat-seq** |

The three papers represent three distinct philosophies: Stutheit-Zhao focuses on biological signal scoring without batch correction; Wilson focuses on absolute technical normalization; Zeng focuses on cross-study harmonization of relative counts. None of the three approaches has been validated against the others in a controlled comparison using the same samples and outcomes, leaving the field without a principled basis for choosing among them.

---

*Review prepared April 2026*
