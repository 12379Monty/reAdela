# Adela vs GRAIL MCED Comparison - RevealJS Presentation

This directory contains a RevealJS presentation comparing Adela's cfMeDIP-seq and GRAIL's targeted bisulfite sequencing platforms for multi-cancer early detection (MCED).

## Files

- `Adela_vs_GRAIL_MCED.Rmd` - Main R Markdown presentation file
- `custom.css` - Custom styling for enhanced appearance
- `README.md` - This file

## Requirements

To render this presentation, you need:

```r
install.packages("revealjs")
install.packages("rmarkdown")
```

## Rendering the Presentation

### Method 1: RStudio
1. Open `Adela_vs_GRAIL_MCED.Rmd` in RStudio
2. Click the "Knit" button
3. The HTML presentation will open in your browser

### Method 2: R Console
```r
rmarkdown::render("Adela_vs_GRAIL_MCED.Rmd")
```

### Method 3: Command Line
```bash
Rscript -e "rmarkdown::render('Adela_vs_GRAIL_MCED.Rmd')"
```

## Presentation Structure

The presentation includes the following sections:

1. **Executive Summary** - Overview of both platforms
2. **Technology Platforms** - Technical approaches and key differences
3. **Detection Performance** - Clinical validation data and sensitivity/specificity metrics
4. **Algorithm & Machine Learning** - Feature engineering and classification methods
5. **Clinical Status** - Regulatory status and commercial availability
6. **Technical Comparison** - Side-by-side feature comparison
7. **Data Availability** - Public datasets and reproducibility resources
8. **Strategic Implications** - Technology maturity and critical evaluation
9. **Conclusions** - Key takeaways and future directions

## Navigation

- **Arrow keys**: Navigate between slides
- **Space bar**: Advance to next slide
- **Esc**: Overview mode (see all slides)
- **S**: Speaker notes mode
- **F**: Fullscreen mode

## Customization

### Change Theme
Edit the YAML header in the .Rmd file:
```yaml
theme: sky  # Options: beige, blood, moon, night, serif, simple, solarized
```

### Modify Styling
Edit `custom.css` to adjust:
- Table colors and formatting
- Font sizes and spacing
- Blockquote styling
- Column layouts

### Add Speaker Notes
Add notes below slides using:
```markdown
<div class="notes">
Your speaker notes here (visible in presenter mode with 'S' key)
</div>
```

## Content Highlights

**Key Technical Comparisons:**
- Bisulfite conversion vs immunoprecipitation
- Targeted vs genome-wide coverage
- Input DNA requirements (1-10 ng vs up to 75 ng)
- Detection sensitivity and specificity metrics

**Clinical Evidence:**
- Adela: CAMPERR study (5,000+ participants), validated MRD
- GRAIL: PATHFINDER 2 (35,878 participants), NHS-Galleri (140,000 participants)

**Data Resources:**
- GEO datasets (GSE79838)
- EGA controlled access (EGAD50000000652, EGAD00001011312)
- GitHub repositories (MEDIPIPE, cfMeDIP-seq-analysis-pipeline)
- Zenodo archives with ML models and reproducibility data

## Output Format

The rendered presentation will be a standalone HTML file that can be:
- Viewed in any modern web browser
- Shared via email or web hosting
- Presented locally without internet connection (all resources embedded)

## Troubleshooting

**Issue: "package 'revealjs' not found"**
```r
install.packages("revealjs")
```

**Issue: CSS not loading**
- Ensure `custom.css` is in the same directory as the .Rmd file
- Check the YAML header includes: `css: custom.css`

**Issue: Slides too crowded**
- Use the `{.smaller}` class on slide headers
- Break content into additional slides
- Adjust font size in custom.css

## License & Attribution

Content derived from comprehensive technical comparison document (November 2025).

**Data Sources:**
- Academic publications (Nature, Annals of Oncology)
- Company press releases and documentation
- Public datasets (GEO, EGA, Zenodo)

**Disclaimer:** Neither Adela nor GRAIL tests are FDA approved/cleared. Information provided for research and educational purposes. Clinical decisions should be made in consultation with healthcare providers.

## Contact

For questions about the source document or cfMeDIP-seq methodology, refer to:
- Princess Margaret Cancer Centre publications
- Adela Inc. (adelabio.com)
- GRAIL (grail.com)
