# cfMeDIP-seq Evolution - RevealJS Presentation

This directory contains a comprehensive RevealJS presentation tracing the evolution of Adela's cfMeDIP-seq (cell-free methylated DNA immunoprecipitation sequencing) methodology from foundational academic protocol to clinical-grade assay platform.

## Files

- `cfMeDIP_Evolution.Rmd` - Main R Markdown presentation file
- `custom_cfmedip.css` - Custom styling for enhanced appearance
- `README_cfMeDIP.md` - This file

## Requirements

To render this presentation, you need:

```r
install.packages("revealjs")
install.packages("rmarkdown")
install.packages("ggplot2")  # For timeline visualization
install.packages("dplyr")    # For data manipulation
```

## Rendering the Presentation

### Method 1: RStudio
1. Open `cfMeDIP_Evolution.Rmd` in RStudio
2. Click the "Knit" button
3. The HTML presentation will open in your browser

### Method 2: R Console
```r
rmarkdown::render("cfMeDIP_Evolution.Rmd")
```

### Method 3: Command Line
```bash
Rscript -e "rmarkdown::render('cfMeDIP_Evolution.Rmd')"
```

## Presentation Structure

The presentation is organized into **10 major sections** with **40+ slides**:

### 1. Protocol Foundation
- cfMeDIP-seq protocol overview (Nature Protocols 2019)
- Four-stage protocol workflow
- Sources of technical variability
- Critical quality control points

### 2. Variability Mitigation
- Systematic vs random variability
- Protocol improvements timeline (2018-2024)
- Batch effect management
- Spike-in normalization (Wilson et al. 2022)

### 3. MEDIPIPE Pipeline
- Automated computational pipeline
- Quality control modules
- Methylation quantification strategies
- Sample aggregation and normalization

### 4. Classifier Evolution
- Foundational classifier (Shen et al. 2018)
- Feature selection methods comparison
- Application-specific classifiers (RCC, breast cancer, multi-cancer)
- Best practices from comparative studies

### 5. Platform Comparisons
- cfMeDIP-seq vs bisulfite sequencing
- cfMeDIP-seq vs targeted methylation methods
- Use case optimization

### 6. Sample Processing
- Physical sample format transitions
- Typical batch structure
- Randomization strategies
- Critical batch boundaries

### 7. Best Practices
- Protocol execution guidelines
- Computational analysis recommendations
- Quality control metrics
- Future research priorities

### 8. Key Publications
- Foundational papers (Nature, Nat Protoc)
- Disease-specific applications (brain, RCC, breast)
- Computational tools (MEDIPIPE)
- Methodology comparisons

### 9. Public Data Resources
- GEO datasets (GSE79838)
- EGA controlled access data
- Zenodo machine learning models
- GitHub analysis repositories

### 10. Evolution Summary
- Interactive timeline visualization
- Key takeaways
- Future outlook (near-term, mid-term, long-term)

## Navigation

- **Arrow keys**: Navigate between slides
- **Space bar**: Advance to next slide
- **Esc**: Overview mode (see all slides)
- **S**: Speaker notes mode
- **F**: Fullscreen mode

## Key Features

### Visual Elements
- **Timeline graph**: R ggplot2 visualization showing evolution from 2018-2024
- **Comparison tables**: Side-by-side methodology comparisons
- **Color-coded variability**: High (red), Moderate (orange), Low (green)
- **Stage diagrams**: Protocol workflow visualization

### Content Highlights

**Technical Innovations:**
- Filler DNA carrier molecule (key innovation)
- Spike-in normalization (5% → <1% batch variance)
- MEDIPIPE automation pipeline
- Multiple classifier approaches (GLMnet, Random Forest, Logistic Regression)

**Clinical Translation:**
- MRD validation (head & neck, lung cancer)
- MCED development (12 cancers, 50+ planned)
- Performance metrics (92-97% sensitivity)
- Regulatory pathway considerations

**Data Resources:**
- Public datasets with accession numbers
- GitHub repositories for reproducibility
- Machine learning models (Zenodo)
- Controlled access patient data (EGA)

## Customization

### Change Theme
Edit the YAML header in the .Rmd file:
```yaml
theme: sky  # Options: beige, blood, moon, night, serif, simple, solarized
```

### Modify Styling
Edit `custom_cfmedip.css` to adjust:
- Table colors and formatting (blue color scheme)
- Variability color coding (red/orange/green)
- Timeline emphasis styles
- Data badge appearances
- Protocol step highlighting

### Add Speaker Notes
Add notes below slides using:
```markdown
<div class="notes">
Your speaker notes here (visible in presenter mode with 'S' key)
</div>
```

## Content Sources

**Primary Publication References:**

1. **Shen et al. 2018** - *Nature* 563:579-583
   - Foundational cfMeDIP-seq methodology
   - GLMnet classifier approach

2. **Shen et al. 2019** - *Nature Protocols* 14:2749-2780
   - Standardized reproducible protocol
   - 3-4 day timeline workflow

3. **Wilson et al. 2022** - *Cell Reports Methods*
   - Spike-in normalization methodology
   - Batch effect reduction validation

4. **Zeng et al. 2023** - *Bioinformatics* 39:btad423
   - MEDIPIPE automated pipeline
   - End-to-end analysis workflow

5. **Halla-aho & Lähdesmäki 2022** - *BMC Bioinformatics* 23:138
   - Comprehensive classifier comparison
   - Best practices recommendations

**Disease-Specific Applications:**

6. **Nassiri et al. 2020** - *Nature Medicine*
   - Brain tumor detection and discrimination

7. **Nuzzo et al. 2020** - *Nature Medicine* 26:1041-1043
   - Renal cell carcinoma detection

8. **De Pascali et al. 2024** - *J Transl Med* 22:934
   - Breast cancer diagnosis and BRCA1/2 detection

## Technical Details

### Protocol Highlights
- **Input DNA:** 1-10 ng cfDNA (ultra-low input)
- **Timeline:** 3-4 days (standard lab)
- **Key innovation:** Filler DNA carrier
- **Critical step:** Immunoprecipitation (overnight)

### Variability Sources
- **High:** Filler composition, antibody lots, IP efficiency
- **Moderate:** Adapter ligation, PCR cycles, sequencing depth
- **Low:** End-repair, purification, size selection

### Classifier Approaches
- **Primary:** GLMnet (elastic net regularization)
- **Alternative:** Random Forest, Logistic Regression
- **Feature selection:** Fisher's exact test, ISPCA
- **Performance:** 92-97% sensitivity across cancer stages

### Data Resources
- **GEO:** GSE79838 (cell line validation data)
- **EGA:** EGAD50000000652 (healthy controls)
- **EGA:** EGAD00001011312 (INSPIRE immunotherapy study)
- **Zenodo:** Multiple ML model repositories

## Output Format

The rendered presentation will be a standalone HTML file that can be:
- Viewed in any modern web browser
- Shared via email or web hosting
- Presented locally without internet connection
- Exported to PDF (via browser print function)

## Use Cases

This presentation is suitable for:
- **Technical audiences:** Bioinformaticians, molecular biologists
- **Clinical researchers:** Understanding methodology evolution
- **Industry stakeholders:** Platform development insights
- **Regulatory discussions:** Technical foundation for clinical translation
- **Educational settings:** Graduate-level molecular diagnostics courses

## Troubleshooting

**Issue: Slides extend beyond the bottom of the screen**
- This has been fixed in `custom_cfmedip.css` with scrolling enabled
- Each slide now scrolls vertically if content is too long
- Font sizes have been reduced for better fitting
- If issues persist, try:
  - Zooming out in browser (Ctrl/Cmd + minus)
  - Using fullscreen mode (F key)
  - Breaking long slides into multiple slides

**Issue: "package 'ggplot2' not found"**
```r
install.packages("ggplot2")
install.packages("dplyr")
```

**Issue: Timeline graph not rendering**
- Ensure ggplot2 and dplyr are installed
- Check R version (>= 4.0 recommended)

**Issue: CSS not loading**
- Ensure `custom_cfmedip.css` is in the same directory as the .Rmd file
- Check the YAML header includes: `css: custom_cfmedip.css`

**Issue: Slides too crowded**
- Use the `{.smaller}` class on slide headers
- Break content into additional slides
- Adjust font size in custom_cfmedip.css

## Technical Evolution Timeline

```
2018: Foundational methodology (Nature)
  ↓
2019: Standardized protocol (Nat Protoc)
  ↓
2020: Disease applications (Brain, RCC)
  ↓
2022: Spike-in normalization (batch ↓5% → <1%)
  ↓
2023: MEDIPIPE automation
  ↓
2024: Clinical validation (MRD, MCED)
  ↓
Future: FDA approval, population screening
```

## Key Takeaways

**Platform Advantages:**
1. Preserves DNA integrity (no bisulfite)
2. Low input requirements (1-10 ng)
3. Genome-wide discovery capability
4. Platform flexibility (same assay, different classifiers)
5. Validated clinical applications (MRD, MCED)

**Evolution Highlights:**
- Batch variance reduction: 5% → <1%
- Automation: Manual → MEDIPIPE pipeline
- Application: Single cancer → multi-cancer platform
- Status: Research tool → clinical validation

## Future Directions

**Near-term (1-2 years):**
- Complete CAMPERR study (5,000+ participants)
- Expand MRD to additional cancer types
- Regulatory pathway progress

**Mid-term (3-5 years):**
- FDA approval for MRD applications
- MCED clinical utility demonstration
- Integration into care pathways

**Long-term (5+ years):**
- Population screening programs
- Multi-omics integration
- Personalized monitoring

## License & Attribution

Content derived from comprehensive technical review of cfMeDIP-seq evolution (January 2026).

**Data Sources:**
- Peer-reviewed publications (Nature, Nature Protocols, etc.)
- Public datasets (GEO, EGA, Zenodo)
- Open-source software (GitHub repositories)

**Disclaimer:** Information provided for research and educational purposes. Clinical decisions should be made in consultation with healthcare providers.

## Contact

For questions about:
- **Original methodology:** Princess Margaret Cancer Centre publications
- **MEDIPIPE pipeline:** `github.com/pughlab/MEDIPIPE`
- **Clinical applications:** Adela Inc. (adelabio.com)
- **Source document:** https://readela.netlify.app/ai/_cfmedip-seq_protocol_pipeline_classifiers

---

*Last updated: January 2026*
