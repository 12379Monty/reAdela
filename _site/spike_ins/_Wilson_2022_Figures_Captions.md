# Figures and Captions from Wilson et al. 2022
## Sensitive and reproducible cell-free methylome quantification with synthetic spike-in controls
### Cell Reports Methods, September 9, 2022

---

## Graphical Abstract

![Graphical Abstract](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/18a8/9499995/c9cf3f5c916a/fx1.jpg)

**Caption:** Wilson et al. developed synthetic DNA controls for normalizing cell-free methylation DNA immunoprecipitation sequencing (cfMeDIP-seq) assays used in liquid biopsies. These spike-in controls allow absolute quantification of methylated cell-free DNA in picomoles rather than arbitrary read counts. They also adjust for batch effects and reduce systematic noise related to technical biases.

---

## Figure 1 - Experimental Design

![Figure 1](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/18a8/9499995/a78845db6b5e/gr1.jpg)

**Figure 1. Experimental design using synthetic spike-in control DNA**

(A) Technical assessment of the spike-in controls with cfDNA mimic. 
  - (Left) Assessment of technical bias in solely the spike-in controls. 
  - (Right) Optimization of the synthetic DNA amount using sheared HCT116 cfDNA mimic.

(B) Clinical evaluation of acute myeloid leukemia (AML) patient samples with spike-in controls.

---

## Figure 2 - Technical Biases Assessment

![Figure 2](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/18a8/9499995/17d0f592d7cc/gr2.jpg)

**Figure 2. Assessing biases in fragment length, G + C content, and CpG fraction**

(A–C) In (A) input spike-in control DNA without cfMeDIP-seq, (B) output spike-in control DNA, after cfMeDIP-seq, and (C) 0.01 ng spike-in control DNA added to HCT116 replicates.

Blue, methylated fragments; gray, unmethylated fragments. Circle, sample 1; triangle, sample 2. Solid line, mean of the two samples. Columns marked with numerals 1 and 2 represent alternative sets of fragments with identical properties but different sequences.

See also Table S2.

---

## Figure 3 - Read Distribution in Genomic Windows

![Figure 3](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/18a8/9499995/6cd01a2af674/gr3.jpg)

**Figure 3. Two-dimensional histograms of the number of reads found in 300 bp windows**

(A and B) Binned by molar amount and either (A) standard deviation of molar amount or (B) Umap k100 multi-read mappability.

Histograms only include windows that do not overlap with UCSC simple repeats and the ENCODE blacklist, and regions with Umap k100 multi-read mappability scores ≤0.5. Asterisks indicate 11 outlier genomic windows chosen for further examination.

---

## Figure 4 - Correlation with EPIC Array

![Figure 4](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/18a8/9499995/4f9569cfa1e5/gr4.jpg)

**Figure 4. Correlation of two measurements of fragment methylation by cfMeDIP and EPIC array M-value for 300 bp genomic windows**

(A, C, E, and G) Molar amount calculated from HCT116 samples correlated to EPIC array M-values.

(B, D, F, and H) Read counts calculated from the same samples, ignoring the spike-in controls.

- (A and B) 37,714 windows with ≥3 CpG probes represented on the EPIC array.
- (C and D) 7,975 windows with ≥5 CpG probes represented on the EPIC array.
- (E and F) 2,066 windows with ≥7 CpG probes represented on the EPIC array.
- (G and H) 158 windows with ≥10 CpG probes represented on the EPIC array.

Solid black line, linear model of best fit; dashed red line, loess (Cleveland 1979) local regression.

---

## Figure 5 - Model Accuracy Assessment

![Figure 5](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/18a8/9499995/3cccb4a8cd6a/gr5.jpg)

**Figure 5. Mean absolute error between known molar amount and predicted molar amount in test data consisting of held-out spike-ins not used for training**

For each number of spike-in fragments between 6 and 25 inclusive, we 100 times randomly selected that number of spike-ins as training data. We used the remaining spike-ins as test data. Each point shows the mean absolute error over all the test spike-ins for that iteration. The vertical limits of the plot include at least 98/100 iterations in every case. We denote outliers for 6 or 7 training spike-ins with a cross at the top of the plot, labeled with the mean absolute error for that case.

Red line denotes median mean absolute error.

See also Table S6.

---

## Figure 6 - Batch Effect Mitigation

![Figure 6](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/18a8/9499995/49611595b438/gr6.jpg)

**Figure 6. Principal component analyses of cfMeDIP results normalized through four different strategies, and associations with experimental variables**

(Left) Proportion of the variance explained by each principal component. 
(Right) Association between known variables, both technical and clinical, and principal component. Cohen's *d* is an effect size of standardized means between variable. ***p < 0.001.

(A) Raw read counts without any normalization.

(B) Read counts normalized using QSEA.

(C) Data normalized using spike-in controls.

(D) Data normalized using spike-in controls and removing regions in UCSC simple repeats, in the ENCODE blacklist, and with Umap k100 multi-read mappability scores ≤0.5.

See also Tables S3, S4, and S5.

---

## Key Figure Insights

### Technical Performance
- **Figure 2** demonstrates that cfMeDIP-seq achieves 97% methylation specificity with ≤0.01% non-specific binding
- Enrichment shows preference for 160 bp fragments (due to size selection), higher G+C content, and higher CpG fraction

### Quantification Accuracy
- **Figure 4** shows strong correlation between spike-in normalized molar amounts and EPIC array M-values (r = 0.62 to 0.82, depending on CpG density)
- **Figure 5** demonstrates high prediction accuracy with mean absolute error ≤0.002 pmol

### Batch Effect Correction
- **Figure 6** shows that spike-in controls reduce batch-associated variance to ≤1% of total variance
- Most effective when combined with filtering of problematic genomic regions

---

## Citation

Wilson SL, Shen SY, Harmon L, Burgener JM, Triche T Jr, Bratman SV, De Carvalho DD, Hoffman MM. Sensitive and reproducible cell-free methylome quantification with synthetic spike-in controls. Cell Rep Methods. 2022 Sep 9;2(9):100294. doi: 10.1016/j.crmeth.2022.100294. PMID: 36160046; PMCID: PMC9499995.
