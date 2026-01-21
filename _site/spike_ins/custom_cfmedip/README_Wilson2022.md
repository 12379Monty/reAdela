# Wilson et al. (2022) Spike-in Controls - RevealJS Presentation

This directory contains a comprehensive RevealJS presentation on the Wilson et al. (2022) publication describing synthetic spike-in controls for cfMeDIP-seq normalization and absolute quantification.

## Files

- `Wilson2022_SpikeIns.Rmd` - Main R Markdown presentation file
- `custom_wilson.css` - Custom dark theme styling
- `README_Wilson2022.md` - This file

## Publication Details

**Full Citation:** Wilson, S.L., Shen, S.Y., Harmon, L., Burgener, J.M., Triche, T., Bratman, S.V., and De Carvalho, D.D. (2022). Sensitive and reproducible cell-free methylome quantification with synthetic spike-in controls. *Cell Reports Methods* 2(9):100294.

**DOI:** 10.1016/j.crmeth.2022.100294  
**PMID:** 36160046  
**PMCID:** PMC9499995

## Requirements

To render this presentation, you need:

```r
install.packages("revealjs")
install.packages("rmarkdown")
install.packages("ggplot2")  # For visualizations
install.packages("dplyr")    # For data manipulation
install.packages("tidyr")    # For data reshaping
```

## Rendering the Presentation

### Method 1: RStudio
1. Open `Wilson2022_SpikeIns.Rmd` in RStudio
2. Click the "Knit" button
3. The HTML presentation will open in your browser

### Method 2: R Console
```r
rmarkdown::render("Wilson2022_SpikeIns.Rmd")
```

### Method 3: Command Line
```bash
Rscript -e "rmarkdown::render('Wilson2022_SpikeIns.Rmd')"
```

## Presentation Structure

The presentation contains **16 major sections** with **60+ slides**:

### 1. Executive Summary
- The problem: No absolute quantification, technical biases, batch effects
- The solution: 52 synthetic spike-in controls
- Key achievement: Batch variance reduced from 5% to <1%

### 2. The Problem (Pre-2022)
- Limitations of cfMeDIP-seq before spike-ins
- Inadequacy of prior approaches (Arabidopsis, input sequencing, QSEA)

### 3. Spike-in Design
- Factorial design: 3×3×3 = 27 conditions
- Parameters: Fragment length (80, 160, 320 bp), GC content (35%, 50%, 65%), CpG fraction (1/80, 1/40, 1/20 bp)
- Dual sequences: Methylated + unmethylated for each
- Design constraints: No human alignment, no secondary structure
- Visual coverage space diagram

### 4. Technical Validation
- **Experiment 1:** Spike-ins alone (9.99 ng pool)
  - Methylation specificity: 97%
  - Non-specific binding: 3%
  - Bias confirmation
  
- **Experiment 2:** HCT116 optimization
  - Optimal amount: 0.01 ng spike-ins per 10 ng sample
  - Uses <1% of sequencing reads
  - Alternative sequence reproducibility
  
- **Experiment 3:** Multi-lab batch effect assessment
  - 5 AML patients × 3 labs = 15 technical replicates
  - Intentional variations across 7 protocol parameters
  - Spike-ins corrected systematic errors

### 5. GLM Quantification Model
- Mathematical framework (Gaussian GLM equation)
- Window-based quantification (300 bp bins)
- Overlap adjustment formula
- Model performance: R² = 0.93, MAE ≤0.002 pmol

### 6. Biological Validation
- Correlation with EPIC array (HCT116)
- 4,446,375 genomic windows analyzed
- Correlation improves with more array probes: r = 0.82 (≥5 probes)
- Validation against orthogonal platform

### 7. Batch Effect Mitigation
- Four normalization approaches compared
- PCA analysis results
- PC1 variance: Biological signals dominate with spike-in + filtering
- QSEA introduced spurious correlation with filler DNA
- Batch variance reduction: 5% → <1% (interactive chart)
- Repetitive element findings (71% RepeatMasker, 43% Alu)

### 8. Quality Control
- Methylation specificity metric and thresholds
- Recommended filtering pipeline (4 categories)
- High molar amount regions (cancer-relevant genes: MIER2, HCN2, RNF126)
- QC decision rules

### 9. The spiky Package
- Bioconductor package overview
- Core functions: `methylation_specificity()`, `model_glm_pmol()`, prediction pipeline
- Workflow integration with code examples
- Compatibility with MEDIPIPE

### 10. Data & Resources
- GEO accession GSE166259 (public)
- EGA accession EGAS00001005069 (controlled access)
- Zenodo deposits (code and package)
- Table S1: Complete spike-in sequences
- Patent PCT/CA2020/051507 licensed to Adela, Inc.

### 11. Key Innovations
- Methodological advances (4 key innovations)
- Impact metrics visualization
- First synthetic spike-ins for enrichment methylation methods

### 12. Limitations & Future
- Acknowledged limitations (4 categories)
- Future opportunities (technical, applications, clinical)
- Platform extensions and clinical translation

### 13. Impact on Field
- Immediate applications (clinical, research, technology)
- Adoption by Adela (commercial integration)
- Broader implications for enrichment methods and cfDNA field

### 14. Cost-Benefit Analysis
- Implementation costs (one-time and per-sample)
- Return on investment metrics
- Negligible cost vs. benefits

### 15. Critical Assessment
- Strengths (experimental design, statistical rigor, open science)
- Areas for development (technical limitations, validation needs)

### 16. Connection to Adela vs GRAIL
- Platform differentiator analysis
- Quantification comparison table
- Integration with MEDIPIPE pipeline
- Key technical advantage

### Conclusions
- Key achievements (4 major points)
- Technical excellence highlights
- Future outlook (near/mid/long-term)

## Key Highlights

### Scientific Achievement

**Problem Solved:**
- **Arbitrary quantification** → Absolute picomolar units
- **Technical biases** → GLM-based correction
- **Batch effects (~5%)** → Reduced to <1% variance

**Innovation:**
- First synthetic spike-in controls for MeDIP-seq
- Factorial design covering fragment property space
- Open-source implementation (spiky package)

**Validation:**
- 3 independent experiments
- Multi-lab reproducibility (3 labs)
- Correlation with EPIC array (r = 0.82)
- 10/10 technical replicates improved

### Commercial Impact

**Patent & Licensing:**
- PCT/CA2020/051507
- Licensed to Adela, Inc.
- Foundation for MCED and MRD platforms

**Adoption:**
- Integrated into MEDIPIPE pipeline
- Used in Adela's clinical validation studies
- Key differentiator vs. GRAIL platform

### Technical Details

**Spike-in Design:**
- 27 factorial combinations (3×3×3)
- 52 total spike-ins (26 methylated, 26 unmethylated)
- Covers: 80-320 bp length, 35-65% GC, 1/80-1/20 bp CpG

**Performance Metrics:**
- Methylation specificity: 97% (target ≥98%)
- Non-specific binding: ≤0.01%
- Model R²: 0.93
- Mean absolute error: ≤0.002 pmol
- Array correlation: r = 0.82 (≥5 probes)

**Quality Control:**
- Pass: ≥0.98 methylation specificity
- Warning: 0.95-0.98
- Fail: <0.95
- Filtering removes: Simple repeats, ENCODE blacklist, low mappability, high variance

## Visualizations

The presentation includes several interactive ggplot2 visualizations:

1. **Spike-in Coverage Space** (Factorial Design Grid)
   - 3×3 matrix showing all parameter combinations
   - Color-coded by synthesis status
   - CpG fraction labels

2. **Batch Variance Reduction** (Bar Chart)
   - Compares 4 normalization methods
   - Shows 5% → <1% improvement
   - Target line at 1%

3. **Impact Metrics** (Grouped Bar Chart)
   - Before/After comparison
   - 4 key metrics: Batch variance, Array correlation, Methylation specificity, Non-specific binding
   - Target values shown as diamonds

## Customization

### Change Theme
Edit the YAML header in the .Rmd file:
```yaml
theme: moon  # Options: beige, blood, moon, night, serif, simple, solarized, sky
```

**Note:** The custom CSS is optimized for the **moon** theme (dark background).

### Modify Styling
Edit `custom_wilson.css` to adjust:
- Dark theme color palette (blues, greens, oranges, reds)
- Table styling for dark backgrounds
- QC threshold color coding (pass/warning/fail)
- Badge appearances (excellent/good/warning)
- Code block backgrounds
- Alert boxes

### Add Speaker Notes
Add notes below slides using:
```markdown
<div class="notes">
Your speaker notes here (visible in presenter mode with 'S' key)
</div>
```

## Navigation

- **Arrow keys**: Navigate between slides
- **Space bar**: Advance to next slide
- **Esc**: Overview mode (see all slides)
- **S**: Speaker notes mode
- **F**: Fullscreen mode
- **B**: Blackout screen

## Content Deep Dive

### Factorial Design Explained

**27 Combinations:**
```
Length (3) × GC (3) × CpG (3) = 27

Length:  80 bp,  160 bp,  320 bp
GC:      35%,    50%,     65%
CpG:     1/80,   1/40,    1/20 bp
```

**Example Combination:**
- 160 bp, 50% GC, 1/40 bp CpG
- Generates 2 sequences: 1 methylated, 1 unmethylated
- Total: 160 bp × 50% GC = 80 bp GC content
- CpG sites: 160 bp / 40 bp/CpG = 4 CpG sites

### GLM Model Equation

**Full Mathematical Model:**

$$\eta = \beta_0 + \beta_{reads} \cdot x_{reads} + \beta_{len} \cdot x_{len} + \beta_{GC} \cdot x_{GC} + \beta_{CpG} \cdot (x_{CpG})^{1/3}$$

**Features:**
- $x_{reads}$: Deduplicated read count (UMI-based)
- $x_{len}$: Fragment length in base pairs
- $x_{GC}$: GC content (0-1 proportion)
- $x_{CpG}$: CpG fraction (cube root transformed)

**Output:**
- $\eta$: Absolute molar amount in picomoles

**Overlap Adjustment:**

$$\eta' = \left(\frac{\ell}{\ell^*}\right) \times \eta$$

Where $\ell$ = overlap length, $\ell^* = 300$ bp (window size)

### Multi-Lab Experiment Details

**5 AML Patients × 3 Labs = 15 Samples**

**Intentional Variations:**

| Parameter | Variation | Purpose |
|-----------|-----------|---------|
| Adapter | xGen Stubby vs Custom | UMI design |
| Ligation | 4°C 16h vs 20°C 2h | Temperature/time |
| Antibody | Lot RD004 vs RD001 | Lot-to-lot variation |
| Filler DNA | 50/50 meth/unmeth vs 100% unmeth | **Protocol violation** |
| PCR Cycles | 11 vs 13 vs 13-15 | Amplification |
| Sequencer | NovaSeq vs NextSeq | Platform |
| Read Depth | 60M vs 100M vs 85M | Coverage |

**Result:** Spike-ins successfully corrected all variations, including the protocol violation in Lab C.

### PCA Analysis Results

**Principal Component 1 (80-83% variance):**

| Method | PC1 Associates With | Cohen's d |
|--------|---------------------|-----------|
| Raw counts | Non-significant technical | Small |
| QSEA | **Filler DNA (artifact)** | Large |
| Spike-in only | Non-significant technical | Small |
| Spike-in + filter | **Biological (sample, sex)** | Large |

**Interpretation:**
- QSEA normalization introduced spurious correlation
- Spike-in + filtering correctly identified biological signals

### Data Resources Summary

**Public Access:**
- GEO GSE166259: HCT116 cell line data
- Zenodo 10.5281/zenodo.4683791: Analysis code
- Zenodo 10.18129/B9.bioc.spiky: spiky package
- Table S1: All 52 spike-in sequences

**Controlled Access:**
- EGA EGAS00001005069: AML patient raw data
- Contact: University Health Network Genomics Data Access Committee

**Software:**
- Bioconductor: `https://bioconductor.org/packages/spiky`
- GitHub: `https://github.com/trichelab/spiky`
- Analysis repo: `https://github.com/hoffmangroup/2020spikein`

## Use Cases

This presentation is suitable for:

**Technical Audiences:**
- Bioinformaticians implementing normalization pipelines
- Molecular biologists optimizing cfMeDIP-seq protocols
- Computational biologists developing methylation analysis methods

**Research Scientists:**
- Understanding spike-in control methodology
- Planning multi-center cfDNA studies
- Evaluating normalization strategies

**Industry Stakeholders:**
- Assay development for clinical diagnostics
- Regulatory submission preparation
- Technology comparison and evaluation

**Clinical Researchers:**
- Understanding Adela platform technology
- Evaluating MCED/MRD assay capabilities
- Multi-center trial design

**Educational Settings:**
- Graduate-level epigenetics courses
- Liquid biopsy methodology seminars
- Bioinformatics training programs

## Troubleshooting

**Issue: Slides extend beyond the bottom of the screen**
- This has been fixed in `custom_wilson.css` with scrolling enabled
- Each slide now scrolls vertically if content is too long
- Font sizes have been reduced for better fitting
- If issues persist, try:
  - Zooming out in browser (Ctrl/Cmd + minus)
  - Using fullscreen mode (F key)
  - Breaking long slides into multiple slides

**Issue: "package 'tidyr' not found"**
```r
install.packages("tidyr")
```

**Issue: Visualizations not rendering**
- Ensure ggplot2, dplyr, and tidyr are installed
- Check R version (>= 4.0 recommended)
- Verify data frame structure in code chunks

**Issue: CSS not loading / Dark theme not applied**
- Ensure `custom_wilson.css` is in the same directory
- Check YAML header includes: `css: custom_wilson.css`
- Verify moon theme is selected: `theme: moon`

**Issue: Math equations not rendering**
- Ensure MathJax is enabled (RevealJS includes by default)
- Check LaTeX syntax in equations
- Try different browser if issues persist

**Issue: Slides too crowded**
- Use `{.smaller}` class on slide headers for extra compact slides
- Break complex content into additional slides
- Adjust font sizes in custom_wilson.css (look for font-size properties)
- Modify base font size: `.reveal { font-size: 28px !important; }` (smaller = 28px, default = 32px)

**Additional Customization:**

To further reduce content size, edit `custom_wilson.css`:

```css
/* Make everything smaller */
.reveal {
  font-size: 28px !important;  /* Reduce from 32px */
}

/* Make tables even smaller */
.reveal table {
  font-size: 0.5em !important;  /* Reduce from 0.55em */
}

/* Tighter line spacing */
.reveal p, .reveal li {
  line-height: 1.2 !important;  /* Reduce from 1.3 */
}
```

## Technical Implementation Notes

**Spike-in Synthesis:**
- 80 bp & 160 bp: Ultramer DNA Oligonucleotides (4 nmol)
- 320 bp: gBlocks Gene Fragments (250 ng)
- Methylation: M.SssI CpG methyltransferase (37°C, 30 min)
- Verification: Restriction digest (HpyCH4IV, HpaII, or AfeI)

**Optimal Spike-in Amount:**
- 0.01 ng per 10 ng biological sample
- <1% of sequencing reads
- >650,000 spike-in reads for GLM training
- >99% reads for biological analysis

**Computational Requirements:**
- R >= 3.6 (4.0+ recommended)
- Bioconductor packages: spiky, GenomicRanges
- Alignment: Bowtie2
- Deduplication: UMI-tools
- RAM: 16 GB minimum, 32 GB recommended

## Key Equations Reference

**1. GLM Quantification:**
```
η = β₀ + β_reads·x_reads + β_len·x_len + β_GC·x_GC + β_CpG·(x_CpG)^(1/3)
```

**2. Window Overlap Adjustment:**
```
η' = (ℓ / ℓ*) × η
```

**3. Methylation Specificity:**
```
MS = N_methylated / N_total
```

**4. Quality Thresholds:**
```
Pass:    MS ≥ 0.98
Warning: 0.95 ≤ MS < 0.98
Fail:    MS < 0.95
```

## Connection to Other Documentation

**Related Presentations:**
1. `Adela_vs_GRAIL_MCED.Rmd` - Platform comparison
2. `cfMeDIP_Evolution.Rmd` - Protocol and classifier evolution

**Integration Points:**
- Wilson 2022 provides normalization layer for cfMeDIP-seq
- Foundational technology for Adela's platform
- Integrated into MEDIPIPE computational pipeline
- Key differentiator in MCED platform comparison

**Documentation Updates:**
- Add spike-in protocol to standard cfMeDIP-seq workflow
- Include GLM quantification in computational pipelines
- Incorporate QC metrics in analysis recommendations
- Update batch effect mitigation strategies

## Citation Information

**Primary Citation:**
```
Wilson, S.L., Shen, S.Y., Harmon, L., Burgener, J.M., Triche, T., 
Bratman, S.V., and De Carvalho, D.D. (2022). Sensitive and reproducible 
cell-free methylome quantification with synthetic spike-in controls. 
Cell Reports Methods 2(9):100294. doi:10.1016/j.crmeth.2022.100294
```

**BibTeX:**
```bibtex
@article{wilson2022sensitive,
  title={Sensitive and reproducible cell-free methylome quantification with synthetic spike-in controls},
  author={Wilson, Samantha L and Shen, Shu Yi and Harmon, Liane and Burgener, Jenna M and Triche, Timothy and Bratman, Scott V and De Carvalho, Daniel D},
  journal={Cell Reports Methods},
  volume={2},
  number={9},
  pages={100294},
  year={2022},
  publisher={Elsevier},
  doi={10.1016/j.crmeth.2022.100294},
  pmid={36160046},
  pmcid={PMC9499995}
}
```

## License & Disclaimer

**Patent Information:**
- PCT/CA2020/051507: Synthetic spike-in controls for cell-free MeDIP sequencing
- Licensed to Adela, Inc.

**Research Use:**
- Open-source spiky package available for research
- Commercial use may require licensing consideration

**Disclaimer:**
- Information for research and educational purposes
- Patent-protected technology
- Contact Adela, Inc. for commercial applications
- Clinical decisions should be made in consultation with healthcare providers

## Output Format

The rendered presentation will be a standalone HTML file that can be:
- Viewed in any modern web browser
- Shared via email or web hosting
- Presented locally without internet connection
- Exported to PDF via browser print function (Ctrl+P, "Save as PDF")

## Summary

This presentation provides a comprehensive technical overview of the Wilson et al. (2022) spike-in normalization methodology, covering:

✓ Problem identification and motivation  
✓ Factorial design strategy (52 spike-ins)  
✓ Three validation experiments  
✓ GLM quantification framework  
✓ Batch effect mitigation (5% → <1%)  
✓ Biological validation (r = 0.82 array correlation)  
✓ Quality control metrics and thresholds  
✓ spiky R package implementation  
✓ Commercial adoption by Adela  
✓ Platform differentiation vs. GRAIL  
✓ Future directions and opportunities

**Key Achievement:** Enables absolute quantification and reduces batch variance by 5-fold, providing foundation for clinical cfMeDIP-seq assays.

---

*Last updated: January 2026*

*Source document: https://readela.netlify.app/spike_ins/_wilson_2022_spike_in_controls_technical_review*
