---
title: Supplementary notes — DNA methylation dynamics
---


Companion to Section 1 of the introductory series. Four questions on the
*dynamics* of the mark — how it comes off, how it copies, how it is
inherited at birth, and what is happening in adult somatic tissue.

---

## 1. The reverse reaction — from methylated to unmethylated

There is no enzyme that simply pops a methyl group off 5-methylcytosine and returns it to cytosine. Biology took the long route. There are two routes off the methylated state:

### Passive demethylation

Replication-coupled dilution. After every round of DNA replication each previously fully-methylated CpG passes through a hemimethylated intermediate (methyl group on the parental strand only). If maintenance methylation fails — i.e. if DNMT1 is not delivered to the daughter strand — the methylated fraction at that locus halves with the next division:

$$\text{fully methylated} \;\to\; \text{hemimethylated} \;\to\; \text{unmethylated}$$

No enzyme is required. Passive loss is the dominant route during pre-implantation development, where DNMT1 is actively excluded from the maternal nucleus.

### Active demethylation

A multi-step oxidation, excision, and repair cascade. The DNA is not directly de-methylated; instead, the modified base is oxidised and then physically replaced via base excision repair.

Step 1 — TET-mediated oxidation. **TET1, TET2, and TET3** are Fe²⁺ / α-ketoglutarate-dependent dioxygenases that iteratively oxidise 5mC:

$$\text{5mC} \;\xrightarrow{\text{TET}}\; \text{5hmC} \;\xrightarrow{\text{TET}}\; \text{5fC} \;\xrightarrow{\text{TET}}\; \text{5caC}$$

Step 2 — base excision. **TDG** (thymine DNA glycosylase) excises 5fC and 5caC, leaving an abasic site.

Step 3 — repair. Base excision repair fills in unmodified cytosine.

### The asymmetry worth noticing

Methylation is one SAM-dependent enzymatic step. Demethylation requires an entire oxidation + glycosylase + repair pathway. The system is mechanistically biased toward **retaining** the mark — a feature, not a bug, given that the methylome is part of the cell's identity record.

### A relevant aside for MeDIP

Anti-5mC antibodies cross-react with 5hmC at variable, lot-dependent levels. Tissues with active TET signalling — particularly brain, where 5hmC reaches ~0.5–1% of cytosines — produce a MeDIP signal that mixes 5mC and 5hmC. This is the basis for the "5mC vs 5hmC indistinguishability" caveat that recurs through Section 3.

---

## 2. How is the methylation mark replicated when cells divide?

DNA replication is semiconservative, so the methyl mark cannot simply be copied like a base. Instead, replication produces a **hemimethylated intermediate** at every previously methylated CpG, and a dedicated machinery re-symmetrises it.

### The reader–writer relay

The CpG dinucleotide is palindromic on the double helix:

```
5'... C p G ...3'        5'... C p G ...3'
       |                         |
       Me                        Me
3'... G p C ...5'   →    3'... G p C ...5'
       |                         (no Me yet on daughter strand)
       Me
```

After replication, the parental strand still carries 5mC; the daughter strand does not. The recovery machinery:

1. **UHRF1** recognises the hemimethylated CpG via its SRA domain. It also reads **H3K9me3** on adjacent histones via its TTD/PHD domains — so the system is coupled to chromatin state, not just to the DNA.
2. UHRF1 **ubiquitinates histone H3**, which (together with PCNA at the replication fork) recruits **DNMT1**.
3. DNMT1 has roughly a 30–40-fold preference for hemimethylated over unmethylated CpG. It methylates the daughter-strand cytosine and restores the symmetric pattern.

### Fidelity

Per-CpG, per-division fidelity is high but not perfect — typical estimates fall in the range of ~96–99%. Two consequences:

- The somatic methylome is **clonally inherited** through cell divisions in a recognisable form.
- The residual per-division error rate is the molecular basis for **epigenetic clocks**. Methylation drift accumulates with mitotic age in a slow but locus-specific way, which is what Horvath, Hannum, PhenoAge, GrimAge and related clocks read out.

### Coupling to chromatin

Maintenance methylation is not a stand-alone DNA process. UHRF1's dependence on H3K9me3 means that loss of the histone mark can drag DNA methylation down with it, and vice versa. This is one mechanism by which heterochromatin domains stay heterochromatic across divisions, and one reason cancer-associated chromatin disorganisation often reads out as methylome drift even without direct mutation of DNA-methylation enzymes.

---

## 3. Inheritance at birth — what is actually established and how

The methylation landscape of a newborn's somatic cells is **not** directly inherited from the parents' somatic methylomes. It is the product of two prior waves of genome-wide reprogramming, with a small set of regions deliberately exempted.

### Wave 1 — pre-implantation reprogramming

Sperm cfDNA arrives heavily methylated (~85% of CpGs). Egg DNA is less methylated (~40%). Immediately after fertilisation:

- The **paternal genome** is demethylated rapidly and **actively**, predominantly TET3-mediated, before the first cleavage division.
- The **maternal genome** is demethylated slowly and **passively**, through replication-coupled dilution as DNMT1 is excluded from the early embryonic nucleus.

The embryo reaches its 5mC nadir around the **blastocyst** stage. Around implantation, **DNMT3A and DNMT3B** then lay down a fresh methylation pattern *de novo*, establishing the embryonic somatic methylome. Lineage-specific elaboration of this pattern continues through gastrulation and organogenesis.

### Wave 2 — germline reprogramming

A second, separate wave of demethylation occurs in **primordial germ cells (PGCs)** as they migrate to the genital ridge. This wave is more thorough than Wave 1 — crucially, it **erases parental imprints**. Sex-specific imprints are then re-established de novo in the developing germline:

- **Male germline:** DNMT3A and DNMT3L re-establish paternal imprints during fetal life.
- **Female germline:** maternal imprints are largely set **post-natally**, during oocyte growth.

### What survives both waves

A small number of regions retain methylation through reprogramming and are inherited in a parent-of-origin manner:

- **Imprinting control regions (ICRs).** ~150 imprinted genes in humans. Their parent-of-origin methylation survives Wave 1 (the entire point of imprinting).
- **Some retrotransposon families.** Notably IAPs in mouse, which retain methylation across reprogramming and provide a substrate for transgenerational epigenetic stability.

### The snapshot at birth

The newborn's somatic cells carry:

- A Wave-1-established, lineage-elaborated **somatic methylome**.
- A complete set of **parent-of-origin imprints**, surviving from Wave 1.

The female germline of the same newborn is **still in the middle of Wave 2** for the next generation — maternal imprint establishment continues post-natally during oocyte growth. This developmental timing is one reason maternal imprinting disorders are biologically plausible at the population level.

### Why "inheritance" is the wrong word for somatic methylation

In the colloquial sense, somatic methylation is **not** inherited from parents the way DNA sequence is. The parental somatic methylomes are not transmitted. What is transmitted is:

- **DNA sequence**, which constrains where DNMT3A/B place de novo marks.
- **Imprints at ICRs**, in a parent-of-origin manner.
- **A small residue** of incompletely-erased marks at specific repeat families.

The rest of the somatic methylome at birth is **established fresh** in the embryo from sequence and chromatin context, not copied from the parents.

---

## 4. Methylation and demethylation in adult somatic cells

In differentiated adult tissue, the methylome is largely stable — but not static. Both methylation and demethylation machineries are active, and the steady state is a dynamic equilibrium.

### Maintenance dominates quantitatively

Every replication forces every previously methylated CpG through a hemimethylated intermediate that the UHRF1–DNMT1 axis must re-symmetrise. In dividing somatic compartments — intestinal crypts, hematopoietic progenitors, epidermal basal layers — this is the dominant methylation activity by far.

### *De novo* activity persists

DNMT3A and DNMT3B remain expressed in most adult tissues, at lower levels than in development. They contribute to:

- Re-establishment of methylation at sites where maintenance has failed.
- Refinement of methylation patterns during cellular differentiation within stem-cell compartments.
- Some context-specific methylation changes in response to signalling (immune activation, hormonal cues).

DNMT3A is the most frequently mutated epigenetic gene in clonal hematopoiesis and AML — a fact more easily understood once one knows that DNMT3A is operationally active in adult HSCs.

### Active demethylation is running at low rate everywhere

5hmC is detectable at ~0.1% of cytosines in most somatic tissues, and substantially higher in brain (~0.5–1%) and embryonic stem cells. This means TET activity is not a developmental relic — it is part of the adult steady state. The functional readouts:

- Continuous low-level erasure at TF-bound enhancers, allowing context-specific gene expression changes.
- Higher activity in tissues with more dynamic transcriptional programs (brain, immune).
- A reservoir of 5hmC that, depending on tissue, may be a stable mark in its own right rather than purely an intermediate.

### Replicative burden varies enormously by tissue

- **High burden:** intestinal crypts (~ daily turnover), hematopoietic compartments, epidermis. The maintenance machinery is under heavy demand; small per-division errors accumulate.
- **Intermediate:** liver, kidney, glandular epithelia.
- **Negligible:** post-mitotic neurons, cardiomyocytes. Here passive drift is essentially absent, but active TET-driven turnover can still operate.

### Aging produces characteristic, locus-specific drift

This is the molecular basis for **epigenetic clocks**. The contributors:

- Imperfect maintenance accumulating over divisions.
- Declining TET activity at some loci, increasing at others.
- Clonal selection in stem-cell compartments — clones with particular methylation signatures expand with age (clonal hematopoiesis is the best-characterised example).
- Senescence-associated chromatin remodeling, which redistributes both 5mC and H3K9me3.

The aggregate signal — measurable as a methylation-derived "biological age" — is reproducible enough across cohorts that it has become a standard covariate in methylation studies.

### Environmental and lifestyle exposures

Reproducible somatic methylation shifts have been documented in response to:

- **Smoking** — *AHRR* and *F2RL3* hypomethylation are among the most reproducible exposure-associated CpG findings in EWAS.
- **Chronic inflammation** — sustained shifts in immune-cell methylomes.
- **Nutritional and metabolic state** — one-carbon metabolism (folate, B12, methionine, choline, betaine) feeds the SAM pool and shifts the SAM:SAH ratio, the methylation-potential of the cell.
- **Age, sex, BMI, alcohol** — all produce reproducible methylation shifts.

These are not minor effects relative to the cancer-associated signals that liquid biopsy is trying to read out from plasma. They are real, characteristic, and a primary reason that case–control studies must match on demographic and lifestyle covariates — the failure mode flagged in Section 2.

### Summary of adult dynamics

| Process | Machinery | Where active | Approximate scale |
|---|---|---|---|
| Maintenance methylation | UHRF1 → DNMT1 | All replicating cells, post-S phase | Largest flux; every CpG, every division |
| *De novo* methylation | DNMT3A, DNMT3B | Stem-cell compartments, differentiation events | Small but locus-specific |
| Active demethylation | TET1/2/3 → TDG → BER | All tissues; highest in brain, ES, immune | ~0.1–1% of cytosines as 5hmC |
| Passive demethylation | DNMT1 failure during replication | Replicating compartments only | Ongoing drift; substrate for aging clocks |

The adult somatic methylome is therefore best thought of as a **cellular identity record**, copied with high but imperfect fidelity at every division, refined by ongoing low-rate writing and erasure, and slowly drifting in characteristic ways with age and exposure. That stability — combined with that drift — is what makes the methylome simultaneously useful as a tissue-of-origin classifier and dangerous as a cancer biomarker if confounders are not controlled.
