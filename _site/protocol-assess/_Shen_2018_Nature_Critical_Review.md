* Critical Review — Shen et al. 2018 (Nature): Proof-of-Concept for cfMeDIP-seq

<!--
**Citation:** Shen SY, Singhania R, Fehringer G, et al. *Sensitive tumour detection and classification using plasma cell-free DNA methylomes.* Nature. 2018;563(7732):579–583. doi:10.1038/s41586-018-0703-0
-->

---

### 1. Framing Considerations {-}

Before evaluating individual tests, two structural points shape how performance claims should be interpreted:

1. **Platform-bound commercialization.** cfMeDIP-seq was developed at Princess Margaret Cancer Centre (De Carvalho lab) and has been continuously positioned for commercialization — first through academic licensing, later through Adela. All publications by the originating group, including this one, must therefore be read as *showcase studies* rather than independent methods evaluations. The usual platform-development workflow is iterative: assay parameters (fragmentation, bead loading, antibody lot, filler DNA composition, library preparation chemistry, sequencing depth, window size, filter thresholds) are tuned on in-house data until the outputs pass the test criteria the authors themselves define. The reported performance is therefore conditional on many unreported design choices.

2. **Assay versioning across studies.** Because the studies are spread over years and the protocol was under active refinement, reported results may reflect *slightly different versions* of the assay. The 2018 Nature paper, the 2019 Nature Protocols paper (Shen et al.), and later Adela publications do not use an identical wet-lab protocol. This complicates any attempt to aggregate performance across publications or to reproduce a single version end-to-end.

3. **Independence assumptions in simulations.** Every *in silico* dilution or spike-in analysis in this paper (and in related papers) models mixed reads as if CpG loci, fragments, and windows were independently and identically drawn. They are not. Real cfDNA signal is correlated across loci through (i) shared fragmentation pattern biases, (ii) CpG-density-dependent enrichment bias, (iii) antibody lot and batch effects that act across loci in a correlated fashion, (iv) GC content and mappability biases, and (v) tumor-fraction-dependent shifts in fragment size distribution. A simulation in which loci behave independently is obviously wrong. Doing something more realistic is hard, but the reported limits of detection derived from such simulations should be treated as *optimistic upper bounds*, not as lower-bound estimates of real-world sensitivity.

<!--
---

-->

### 2. Test-by-Test Summary {-}

The paper contains roughly eight distinct tests, each with its own criteria and data. Listed in approximate order of appearance:

| # | Test | Primary metric |
|---|------|---------------|
| 1 | Low-input protocol validation | Library complexity, specificity of enrichment |
| 2 | cfMeDIP-seq vs WGBS concordance | Correlation across methylation calls |
| 3 | Technical reproducibility | Replicate correlation |
| 4 | Analytical sensitivity — *in vitro* cell-line dilution | Recovery of expected DMRs at decreasing spike-in |
| 5 | Analytical sensitivity — *in silico* AML read dilution | Classifier signal retention at decreasing tumor read fraction |
| 6 | Binary classification (AML vs healthy) | AUROC, sensitivity, specificity |
| 7 | Multi-cancer classification (7 tumor types + healthy) | Per-class AUROC, overall accuracy |
| 8 | Tissue-of-origin / early-stage detection | Per-stage sensitivity, class confusion |

---

### 3. Detailed Per-Test Review {-}

#### Test 1 — Low-input protocol validation {-}

**Performance criteria.** The assay should produce libraries of acceptable complexity and CpG-enrichment specificity across an input-DNA titration, with emphasis on the low end relevant to plasma cfDNA (low ng range).

**Study data description.** Sonicated genomic DNA from reference sources was titrated across inputs typically reported as 1 ng, 5 ng, 10 ng, and 100 ng. Filler DNA (λ) was added to bring total DNA into the workable range for immunoprecipitation. Enrichment specificity was assessed using methylated/unmethylated spike-in controls included in the library prep.

**Study data performance.** The paper reports that libraries prepared from low input (down to ~1–10 ng) achieve CpG enrichment and complexity metrics comparable to higher-input libraries.

**Necessary assumptions and caveats.**

- Sonicated genomic DNA does *not* recapitulate plasma cfDNA fragment-size distribution (mononucleosomal ~167 bp with characteristic di-/tri-nucleosome peaks) or its cell-type-of-origin distribution. Extrapolating "works at 10 ng of sonicated DNA" to "works at 10 ng of plasma cfDNA" is a conditional claim.

- Filler DNA specification is critical: later work (Wilson et al. 2022 and the spiky framework) makes clear that lambda-only filler (fully unmethylated) biases the downstream methylation calls. The 2018 paper predates this refinement; its input-titration conclusions should be re-evaluated against current knowledge that filler should be 50% methylated / 50% unmethylated.

- Enrichment metrics reported are *relative* (methylated vs unmethylated spike-in). They do not quantify absolute bias across CpG-density strata, which is where most downstream modeling error enters.

---

#### Test 2 — cfMeDIP-seq vs WGBS concordance {-}

**Performance criteria.** cfMeDIP-seq methylation signal should correlate with WGBS-derived methylation on matched samples, justifying its use as a whole-methylome readout.

**Study data description.** Matched DNA samples processed by both cfMeDIP-seq and WGBS (typically sonicated genomic DNA from cell lines; sample numbers small). Concordance assessed by Pearson/Spearman correlation of per-window or per-region methylation signal.

**Study data performance.** High reported correlation at the regional level (numbers in the 0.8+ range for appropriately binned windows).

**Necessary assumptions and caveats.**

- Correlation is *not* a measure of quantitative accuracy. cfMeDIP-seq is an enrichment-based readout producing read counts proportional to methylated fragment abundance (modulated by CpG density, fragment length, GC content, and antibody affinity). WGBS produces per-CpG fractional methylation. A high correlation at the region level says the two readouts rank regions similarly; it does not say cfMeDIP-seq counts can be inverted to beta-values.

- The correlation is computed on cell-line or bulk-tissue DNA, not cfDNA, and so does not address whether cfMeDIP-seq of plasma recapitulates cfDNA methylation as it would be measured by cfWGBS.

- Concordance analyses typically exclude regions with low signal in both assays, which inflates the apparent correlation.

- **Visualization gap:** Correlations are reported without accompanying scatterplots. For count-based vs fractional readouts, scatterplots would reveal whether the relationship is approximately linear, whether outliers drive the correlation, and whether a log transformation (with appropriate handling of zeros) would better capture the underlying relationship. The absence of this basic visualization makes it difficult to assess whether Pearson correlation is the appropriate summary statistic.

---

#### Test 3 — Technical reproducibility {-}

**Performance criteria.** Replicate libraries from the same input should produce highly correlated genome-wide signal.

**Study data description.** Technical replicates of cfMeDIP-seq libraries. Numbers are small; reproducibility assessed by correlation.

**Study data performance.** Reported as high (correlation coefficients near 1 for shared reference inputs).

**Necessary assumptions and caveats.**

- Same-run technical replicates share batch, antibody lot, operator, and reagent lot. They are near-upper-bound estimates of reproducibility. Between-batch, between-lot, and between-site reproducibility — the relevant quantities for a clinical assay — are not reported in this paper and remain a known weakness of cfMeDIP-seq assays generally (documented in later work; see Wilson et al. 2022).

**Assessment method:** 

- Reproducibility is summarized using correlation coefficients, which capture the strength of linear association but provide little insight into the structure of disagreement between replicates. 

- For high-throughput sequencing data, the established best practice — dating to the beginning of the omics era — is to visualize reproducibility using MA plots (mean vs difference plots, also known as Bland-Altman plots or Tukey mean-difference plots). 

- MA plots reveal systematic bias between replicates and, critically, show how reproducibility degrades at low signal levels, where most genomic regions reside. 
- The decrease in reproducibility at the low end becomes obvious in such plots, and interpretable summaries (such as limits of agreement across signal ranges) can be derived. 

By reporting only correlation coefficients, the authors bypass this rich analytical tradition and provide a summary that is impossible to interpret beyond ordering the strength of linear relationships.

---

#### Test 4 — Analytical sensitivity: *in vitro* cell-line dilution {-}

**Performance criteria.** The assay should recover expected cell-line-specific methylation signal from samples in which the target DNA has been serially diluted into a non-target background.

**Study data description.** Sonicated HCT116 (colorectal carcinoma) DNA titrated into sonicated MCF10A (near-normal breast epithelial) DNA at decreasing fractions — typically spanning 100%, 10%, 1%, 0.1%, 0.01%, and 0.001%. Total input in the low-ng range, filler DNA added as in Test 1.

**Study data performance.** The paper reports that HCT116-characteristic differential methylation signal is recoverable at very low spike-in fractions, with claims extending into the 0.001% range.

**Necessary assumptions and caveats.**

- **Model mismatch with plasma.** Sonicated cell-line genomic DNA does not model plasma cfDNA. Real cfDNA from a tumor is shorter, is chromatin-positioned, and has a tumor-fraction-dependent size distribution. Signal in plasma is not simply a linear mixture of two bulk genomic profiles.

- **Independence.** Sensitivity at the fraction level is estimated by looking at DMR-level summary statistics that implicitly treat loci as independent draws. Under realistic correlated noise (batch, antibody, GC bias), the effective sample size is much smaller and detection limits are accordingly worse.

- **Choice of discriminating DMRs.** DMRs used to assess recovery were typically selected on the undiluted (100%) samples, then evaluated on the diluted samples. This is circular: the features are guaranteed to carry signal in the undiluted endpoint and the question reduces to whether enough reads remain at dilution to preserve ranking. The genuine question — whether an *a priori* panel calibrated in one cohort detects tumor DNA at 0.001% in unseen plasma samples — is not answered.

- **Reporting convention.** "Detection at 0.001%" usually means "the collective DMR signal is statistically distinguishable from background in a small number of technical replicates," not "any individual sample at 0.001% is correctly classified."


---

#### Test 5 — Analytical sensitivity: *in silico* AML read dilution {-}

**Performance criteria.** Tumor-specific methylation signal should be recoverable when reads from a tumor-bearing sample are computationally mixed into reads from a healthy sample, across decreasing mixing fractions.

**Study data description.** Reads from an AML patient cfMeDIP-seq library were downsampled and mixed with reads from a healthy-control cfMeDIP-seq library at specified fractions. The mixed read sets were then re-processed through the classifier pipeline and scored.

**Study data performance.** Classifier output (methylation score) declines gradually with decreasing tumor fraction; the paper reports retention of a discriminating signal at very low fractions.

**Necessary assumptions and caveats.**

<p><p/>
- **This is the test most directly affected by the independence assumption.** *In silico* mixing of aligned reads treats the two libraries as if they were generated by independent draws from two fixed multinomial distributions over loci. The actual data-generating process is:
  - correlated CpG-density-dependent enrichment across loci,
  - correlated GC and fragment-length biases,
  - batch and antibody-lot effects that shift a non-random subset of loci together,
  - overdispersion within samples that is *not* reproduced by Poisson-like thinning.

<p><p/>
- A simulation in which loci behave independently is obviously wrong. The reported limit of detection is therefore an optimistic upper bound under a counterfactual in which every source of correlated error has been removed by construction. Doing something more realistic — e.g., mixing whole samples generated by repeated independent library preparations from known-tumor-fraction plasma, or using a correlated negative-binomial generative model fit to replicate data — is hard, but is the only principled way to anchor a detection limit.

- The AML-into-healthy mixing also ignores that AML cfDNA is heavily hematopoietic in origin and therefore shares much of its methylation pattern with healthy leukocyte-derived cfDNA. Detection therefore rides on a small number of AML-specific DMRs that are themselves selected using the very samples being mixed.

- Classifier thresholds and feature panels applied to the mixed samples are those trained on the unmixed endpoint. Any *in silico* "validation" of detection limits here is internal; no held-out validation of the mixing sensitivity at each level is performed.

---

#### Test 6 — Binary classification: AML vs healthy {-}

**Performance criteria.** A classifier trained on cfMeDIP-seq from AML patients and healthy controls should discriminate the two classes with high AUROC on a held-out set.

**Study data description.** Plasma cfMeDIP-seq libraries from a modest number of AML patients and healthy controls. Feature space reduced to ~300 bp genomic windows (or binned regions) overlapping CpG islands, shores, shelves, and FANTOM5 enhancers. Classifier: GLMnet logistic regression with cross-validation.

**Study data performance.** AUROC reported near the ceiling.

**Necessary assumptions and caveats.**

- Cohort composition: AML and healthy-control samples were collected and processed under different logistical conditions (different sites, time periods, preanalytical workflows). Batch is confounded with disease status. Any classifier built on such data will partly learn batch rather than disease. The paper does not fully decouple these.

- Feature selection (the top-*k* variable or differential windows) was in most Shen-lab analyses performed on the whole dataset prior to cross-validation, or inside a cross-validation in which information leakage is not prevented at the DMR-selection step. This inflates AUROC estimates relative to a pipeline in which *all* steps (including window selection) are nested inside each fold.

- Sample sizes (tens of patients per class) are small; AUROC estimates have wide confidence intervals that are typically not reported.

---

#### Test 7 — Multi-cancer classification (seven tumor types + healthy) {-}

**Performance criteria.** A one-vs-rest (or one-vs-healthy) classifier trained on cfMeDIP-seq from multiple cancer types should discriminate each cancer type from healthy and, secondarily, from the other cancer types, at high AUROC on held-out data.

**Study data description.** Plasma cfMeDIP-seq from patients across ~7 cancer types (pancreatic, colorectal, breast, lung, bladder, renal, AML) plus healthy controls. Typical sample sizes per class were on the order of tens. Features: 300 bp windows in CpG islands / shores / shelves / FANTOM5 enhancers. Classifier: GLMnet logistic regression; reported in one-vs-rest configuration with cross-validation.

**Study data performance.** High per-class AUROCs are reported, typically in the 0.9+ range for most cancers vs healthy; lower for pairwise discrimination between cancer types.

**Necessary assumptions and caveats.**

- All the feature-selection-leakage concerns from Test 6 apply and are more severe because the feature space is larger and the per-class sample count is smaller.
- Site-of-collection / batch confounding: samples for different cancer types were drawn from different cohorts and different institutions. The classifier may be learning collection-site signal rather than cancer-type signal. This is a well-documented concern in later cfMeDIP-seq work (Stutheit-Zhao et al. 2024 explicitly motivates site correction).

- One-vs-rest with highly imbalanced class counts gives per-class AUROCs that are easier to inflate than pairwise AUROCs; the latter are less prominently reported.

- Generalization to new patients, new collection sites, and new assay runs is not tested in this paper. Later work has repeatedly shown that cfMeDIP-seq classifiers trained in one batch degrade substantially when applied to an independent batch without normalization (motivating the Wilson 2022 spike-in approach and the Stutheit-Zhao 2024 pipeline).


---

#### Test 8 — Tissue-of-origin and early-stage detection {-}

**Performance criteria.** The multi-cancer classifier should correctly assign tissue of origin, and should retain sensitivity in early-stage disease where tumor DNA shedding is lowest.

**Study data description.** Early-stage (stage I/II) cases were available in limited numbers for some cancer types (notably pancreatic). Tissue-of-origin analysis used the same GLMnet-derived feature space; stage stratification was performed post hoc.

**Study data performance.** The paper reports that early-stage cases are detected, including for pancreatic cancer. Tissue-of-origin accuracy is reported per class.

**Necessary assumptions and caveats.**

- Early-stage subgroup sizes are small; point estimates of per-stage sensitivity carry wide confidence intervals that are not prominently reported.

- Stage stratification is post hoc from a classifier *tuned on all stages*. A classifier genuinely optimized for stage I detection would likely select different features, and its stage I performance is not what is being reported.

- Tissue-of-origin accuracy is reported on the same samples used (implicitly) in feature space construction. There is no independent tissue-of-origin validation cohort.

- The broader claim that cfMeDIP-seq is especially suited to low-shedding tumors (brain, renal) rests on subsequent papers (Nassiri 2020, Nuzzo 2020), not on this one.

---

### 4. Overall Assessment {-}

Shen et al. 2018 is a well-executed proof-of-concept that establishes cfMeDIP-seq as a technically feasible genome-wide methylome assay from low-input cfDNA and demonstrates that GLMnet classifiers built on binned-window count features can discriminate cancer from healthy plasma on in-group data. As a **proof of concept**, it succeeds.

As **evidence for clinical performance**, its claims should be read with the following qualifications:

1. Every reported limit of detection from *in silico* or *in vitro* dilution is an optimistic upper bound because the underlying noise model assumes locus independence that the data do not have.

2. Classifier AUROCs are conditional on the specific cohorts, assay version, feature-selection pipeline, and batch structure used; none of the three were independently validated in this paper.

3. Batch and site effects are not disentangled from disease effects in the experimental design; later work in the same lab explicitly concedes that normalization is necessary to address this, implying that the 2018 numbers do include a batch component.

4. Results across this paper and later cfMeDIP-seq publications reflect evolving versions of the wet-lab and analytical pipeline; aggregated performance across studies should not be interpreted as replicated performance of a single assay.

The appropriate interpretation is: cfMeDIP-seq from 2018 is a promising platform whose clinical potential requires validation under (i) pre-specified assay parameters, (ii) batch-balanced multi-site collection, (iii) nested feature selection, and (iv) prospective held-out cohorts. These conditions are not met in Shen et al. 2018, and claims downstream of that paper should be evaluated against whether they meet them.

---

*Review prepared as a technical appraisal, not a rebuttal. The goal is to identify the necessary assumptions behind reported performance so that downstream use of these results is appropriately conditioned.*
