# COSMIC Mutational Signatures in Early-Onset vs Late-Onset Breast Cancer

[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.17+-1a9641?style=flat-square)](https://bioconductor.org/)
[![TCGA](https://img.shields.io/badge/Data-TCGA--BRCA-blue?style=flat-square)](https://portal.gdc.cancer.gov/)
[![COSMIC](https://img.shields.io/badge/Signatures-COSMIC_v3.3.1-orange?style=flat-square)](https://cancer.sanger.ac.uk/signatures/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)]()
[![GDC Client](https://img.shields.io/badge/GDC_Client-v1.6.1-lightgrey?style=flat-square)](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20WSL2-informational?style=flat-square&logo=linux&logoColor=white)]()

---

## Overview

This project investigates differences in **somatic mutational signatures** between **early-onset** (≤45 years) and **late-onset** (≥55 years) breast cancer patients using whole-exome sequencing (WXS) data from **The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA)** cohort.

Mutational signatures reflect the underlying biological processes — DNA repair defects, APOBEC activity, ageing-related damage — that drive tumour evolution. By comparing signature exposures across age groups, this analysis aims to uncover **age-specific mutagenic mechanisms** in breast cancer.

---

## Scientific Background

| Signature | Biological Process | Relevance to BRCA |
|-----------|-------------------|-------------------|
| **SBS3**  | Homologous recombination deficiency (HRD) | BRCA1/2 dysfunction |
| **SBS2**  | APOBEC cytidine deaminase activity | Elevated in younger patients |
| **SBS13** | APOBEC cytidine deaminase activity | Co-occurs with SBS2 |
| **SBS5**  | Clock-like, age-associated damage | Higher in older patients |
| **SBS1**  | Spontaneous deamination of 5-methylcytosine | Age-correlated |

Early-onset breast cancer is often associated with germline BRCA1/2 mutations and HRD (SBS3), while late-onset disease tends to show more clock-like signatures (SBS5, SBS1) due to accumulated replication errors over time.

---

## Repository Structure

```
COSMIC-Signatures-EarlyOnset-BRCA/
│
├── scripts/
│   ├── signature_analysis.R       # Main analysis pipeline
│   ├── plot_signatures.R          # Comprehensive visualisation script
│   ├── fix_prepare.R              # MAF file ingestion utility
│   └── fix_clinical.R             # Clinical data fetching utility
│
├── results/
│   ├── plots/
│   │   ├── plot1_mutation_spectrum.jpg    # 96-channel trinucleotide profiles
│   │   ├── plot2_stacked_bar.jpg         # Per-sample stacked bar chart
│   │   ├── plot3_heatmap.jpg             # Signature exposure heatmap
│   │   ├── plot4_boxplots.jpg            # Signature boxplots by group
│   │   ├── plot5_volcano.jpg             # Volcano plot
│   │   └── plot6_cosine_similarity.jpg   # Reconstruction quality
│   └── data/
│       ├── maf_combined.rds              # Combined MAF data (992 samples)
│       ├── clinical.rds                  # Clinical metadata
│       └── signature_stats_results.csv   # Statistical comparison results
│
├── gdc_tool/
│   ├── gdc-client                        # GDC Data Transfer Tool binary
│   ├── gdc_manifest.txt                  # Download manifest (992 files)
│   └── gdc_client_configuration.dtt      # GDC client config
│
├── .gitignore
└── README.md
```

---

## Data

| Property | Details |
|----------|---------|
| **Project** | TCGA-BRCA |
| **Data type** | Masked Somatic Mutation (WXS) |
| **Workflow** | Aliquot Ensemble Somatic Variant Merging and Masking |
| **Genome build** | GRCh38 / hg38 |
| **Total samples** | 992 |
| **Early-onset (≤45)** | Extracted from GDC clinical API |
| **Late-onset (≥55)** | Extracted from GDC clinical API |
| **COSMIC version** | v3.3.1 SBS signatures |

---

## Dependencies

### R Packages

```r
# Bioconductor
BiocManager::install(c(
    "TCGAbiolinks",
    "maftools",
    "MutationalPatterns",
    "BSgenome.Hsapiens.UCSC.hg38",
    "GenomicRanges",
    "GenomeInfoDb",
    "Biostrings"
))

# CRAN
install.packages(c(
    "ggplot2",
    "dplyr",
    "tidyr",
    "pheatmap",
    "RColorBrewer",
    "gridExtra",
    "data.table",
    "jsonlite"
))
```

### System Tools

- [GDC Data Transfer Tool v1.6.1](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) (Ubuntu x64)
- R ≥ 4.2
- Bioconductor ≥ 3.17

---

## Reproducing the Analysis

### 1. Clone the repository

```bash
git clone https://github.com/g-Poulami/COSMIC-Signatures-EarlyOnset-BRCA.git
cd COSMIC-Signatures-EarlyOnset-BRCA
```

### 2. Download the data

```bash
# Download 992 TCGA-BRCA MAF files using the provided manifest
./gdc_tool/gdc-client download -m gdc_tool/gdc_manifest.txt -n 1
```

> **Note:** The `-n 1` flag forces single-threaded downloading, which avoids a known Python multiprocessing pickling bug in the GDC client on WSL/Linux environments.

### 3. Prepare the MAF data

```bash
Rscript scripts/fix_prepare.R
```

This reads all 992 `.maf.gz` files, handles type inconsistencies across files, and saves a combined `results/data/maf_combined.rds`.

### 4. Fetch clinical data

```bash
Rscript scripts/fix_clinical.R
```

Fetches age-at-diagnosis data directly from the GDC API and saves `results/data/clinical.rds`.

### 5. Run the main analysis

```bash
Rscript scripts/signature_analysis.R
```

### 6. Generate all plots

```bash
Rscript scripts/plot_signatures.R
```

---

## Outputs & Interpretation

### `signature_stats_results.csv`

Statistical comparison of COSMIC signature exposures between early-onset and late-onset groups.

| Column | Description |
|--------|-------------|
| `Signature` | COSMIC SBS signature identifier |
| `Early_median` | Median absolute exposure in early-onset patients |
| `Late_median` | Median absolute exposure in late-onset patients |
| `Log2FC` | Log2 fold change (Early / Late) |
| `W_statistic` | Mann-Whitney U test statistic |
| `P_value` | Raw p-value |
| `P_adj` | Benjamini-Hochberg FDR-corrected p-value |

Signatures with `P_adj < 0.05` and `|Log2FC| > 0.5` are considered statistically and biologically significant.

---

### Plot 1 — 96-Channel Trinucleotide Mutation Spectrum

![96-Channel Trinucleotide Mutation Spectrum by Onset Group](results/plots/plot1_mutation_spectrum.jpg)

**What is shown:** The relative frequency of all 96 single base substitution (SBS) mutation types across six substitution classes (C>A in blue, C>G in black, C>T in red, T>A in grey, T>C in green, T>G in pink), each broken into 16 trinucleotide contexts, averaged across patients in the early-onset (top row) and late-onset (bottom row) groups.

**Detailed observations from this plot:**

- **C>G transversions are the most visually striking feature of the early-onset spectrum.** Several contexts within the C>G panel (most notably T[C>G]T and a flanking context) show very tall black bars in early-onset that are substantially reduced in late-onset. This pattern is consistent with an active APOBEC-like or oxidative mutagenic process operating more prominently in younger patients.
- **C>T transitions dominate both groups**, as expected for breast cancer. The tallest C>T bar across both groups corresponds to the T[C>T]A or similar TpCpN context, a hallmark of APOBEC activity (SBS2/SBS13). This peak appears modestly higher in early-onset relative to late-onset when normalised, suggesting a marginally greater APOBEC contribution in the younger cohort.
- **C>A transversions (blue) show a notable asymmetry.** In the late-onset panel, a single tall blue bar emerges at the rightmost C>A context (likely T[C>A]T), which is absent or minimal in early-onset. This context can reflect SBS18 or oxidative damage accumulating with age.
- **T>C, T>A, and T>G substitutions are uniformly low in both groups**, with no dramatic differences visible between panels, suggesting these processes (e.g., SBS5 clock-like damage, which distributes broadly across T>C contexts) do not create strongly context-specific peaks detectable at this resolution.
- **Overall spectrum similarity:** The two profiles are broadly similar in shape — both are dominated by C>T and C>G — indicating that the major mutational processes are shared between age groups. Differences are quantitative rather than qualitative, which is consistent with the statistical results showing limited significance after FDR correction (see Plot 5).

---

### Plot 2 — Stacked Bar Chart of Signature Contributions

![Relative COSMIC Signature Contributions per Sample](results/plots/plot2_stacked_bar.jpg)

**What is shown:** The relative contribution of each COSMIC SBS signature to every individual sample (one vertical bar = one patient), split into early-onset (left panel) and late-onset (right panel). Each colour represents a distinct COSMIC signature; only signatures with a mean relative exposure >1% are displayed. The y-axis represents relative exposure (values above 1.0 indicate that the refitting algorithm assigned more total weight than the observed mutation count — an artefact of the refitting procedure in low-mutation-burden samples).

**Detailed observations from this plot:**

- **Extreme inter-patient heterogeneity is the dominant feature in both groups.** No single signature uniformly dominates; instead, the colour composition of each bar varies dramatically from patient to patient. This reflects the well-known biological heterogeneity of breast cancer and is consistent with the diverse molecular subtypes (ER+, HER2+, TNBC) pooled together in TCGA-BRCA.
- **Several outlier samples exceed a relative exposure of 1.0–2.3.** These samples, visible as very tall bars in both panels (particularly in late-onset where one sample reaches ~2.3), likely have very low total mutation counts, causing the NNLS/refitting algorithm to over-assign signature weights. These are statistical artefacts and should be treated with caution in downstream interpretation.
- **SBS1 (teal/cyan, clock-like CpG deamination) contributes a consistent baseline across most samples in both groups,** visible as a common colour at the base of many bars. This is expected given that SBS1 accumulates with cell divisions throughout life and is nearly universal in cancer.
- **SBS2 (salmon/pink-red, APOBEC) and SBS13 (coral-red, APOBEC) are visible as prominent coloured segments in a subset of samples in both groups,** not exclusively in early-onset. This suggests APOBEC-driven mutagenesis is heterogeneous and sample-specific rather than age-stratified at the group level.
- **SBS5 (yellow-green, clock-like replication) appears to contribute a substantial yellow component at the base of many bars in both groups,** which is consistent with its role as a ubiquitous background process. It does not appear markedly higher in late-onset as a visual impression.
- **Many minor signatures (SBS21, SBS24, SBS29, SBS39, SBS42, SBS86, SBS87) contribute small fragments** to individual samples and likely reflect sparse, sample-specific noise from the refitting, particularly in samples with low mutation burden.
- **The early-onset panel appears slightly more compact** (fewer extreme outlier bars), while late-onset has more tall outlier spikes, likely because of greater sample size and therefore more extreme low-burden samples in the larger late-onset cohort.

---

### Plot 3 — Signature Exposure Heatmap

![Signature Exposure Heatmap (Early vs Late Onset)](results/plots/plot3_heatmap.jpg)

**What is shown:** A hierarchically clustered heatmap where each column is a patient sample and each row is a COSMIC SBS signature. Colour intensity (white → orange → dark red) represents relative exposure magnitude. The annotation bar at the top of the heatmap colour-codes samples by group: red = early-onset, blue = late-onset. Row and column dendrograms reflect unsupervised hierarchical clustering.

**Detailed observations from this plot:**

- **The column dendrogram does not separate early-onset from late-onset samples into distinct clusters.** The group annotation bar at the top shows a thoroughly interleaved distribution of red and blue samples across the entire dendrogram. This is a critical finding: the global signature profile of a tumour does not reliably predict which age group the patient belongs to. Mutational landscape is not a strong discriminator of onset age at the whole-signature-profile level.
- **SBS2 and SBS13 (the two APOBEC signatures) cluster together at the top of the row dendrogram and show the highest intensity values in a subset of samples.** These high-exposure samples (visible as deep orange-to-red columns in the SBS2/SBS13 rows) are scattered across both early- and late-onset groups, reinforcing that APOBEC activity is a sample-specific rather than age-stratified phenomenon in this cohort.
- **SBS87 clusters immediately below SBS2/SBS13,** suggesting that when APOBEC signatures are high, SBS87 also tends to be elevated in the same samples — possibly due to misattribution during the refitting process, or a genuine co-occurring process.
- **SBS39, SBS7a, SBS10b, SBS42, SBS86, SBS6, SBS24, SBS50, SBS16, SBS54, SBS31, SBS29, SBS8, SBS22, SBS30, SBS19, SBS21, SBS7b, SBS59, SBS32, SBS58** all show predominantly pale/white colouration across almost all samples, indicating very low or near-zero relative exposures. Their presence in the refitting output likely reflects the statistical limits of signature attribution in WXS data with moderate mutation burden.
- **SBS1 and SBS15 cluster at the bottom of the row dendrogram.** SBS1 shows a notable high-exposure cell in one or two early-onset samples (visible as a red square at the bottom-left of the heatmap), and some elevated exposure spread across several late-onset samples, consistent with its age-correlated nature — though the pattern is not uniform.
- **The heatmap confirms that most samples have a signature profile characterised by a small number of active signatures and many near-zero contributions,** which is typical for WXS data where mutation counts per sample are insufficient to reliably decompose more than 2–3 signatures simultaneously.

---

### Plot 4 — Signature Boxplots by Group

![Signature Exposures: Early vs Late Onset BRCA](results/plots/plot4_boxplots.jpg)

**What is shown:** Boxplots comparing the relative exposure of COSMIC signatures that reached statistical significance (FDR < 0.05) between early-onset (red) and late-onset (blue) groups. Individual patient data points are overlaid as jitter. Eight signatures are shown: SBS1, SBS2, SBS30, SBS34, SBS38, SBS54, SBS60, and SBS88.

**Detailed observations from this plot:**

- **SBS1 (clock-like CpG deamination):** Both groups have similar medians (approximately 0.08–0.10), with the interquartile ranges substantially overlapping. The late-onset group appears to have a marginally wider spread and slightly more patients at higher exposures, consistent with the known age-correlated accumulation of SBS1. However, there are outliers with SBS1 relative exposures exceeding 0.75–0.90 in both groups, indicating that SBS1 dominance is sample-specific and not age-exclusive.
- **SBS2 (APOBEC — C>T at TpCpN):** The early-onset median is modestly higher than late-onset, with the early-onset IQR sitting slightly above zero while the late-onset box spans a similar low range. However, both groups show extensive right-tail outliers reaching 0.3–0.49, and the overall difference is small. This is directionally consistent with the hypothesis of greater APOBEC activity in early-onset, but the magnitude is modest.
- **SBS30 (base excision repair deficiency):** Median exposure is effectively zero in both groups. Both distributions are heavily right-skewed with a few outlier samples showing exposures up to 0.2–0.4. No meaningful difference is visually apparent between groups.
- **SBS34 (unknown aetiology):** Both groups show medians at or near zero, with most samples having essentially no exposure. A handful of late-onset outliers reach up to 0.075–0.10. This signature is unlikely to be biologically meaningful at the cohort level.
- **SBS38 (damage from UV light or platinum chemotherapy — context uncertain in breast cancer):** Medians are essentially zero in both groups, with sparse outliers up to 0.30 in late-onset and 0.05 in early-onset. The signature is not biologically expected to be prominent in breast cancer.
- **SBS54 (artefact or unknown):** Late-onset shows a visibly higher median and IQR compared to early-onset, making this the most visually distinct comparison in the panel. However, the absolute exposure values are small, and SBS54 is not a well-characterised biologically meaningful signature in breast cancer — this likely reflects a technical or low-burden refitting artefact.
- **SBS60 and SBS88:** Both have near-zero medians in both groups with sparse high outliers. These are poorly characterised signatures with uncertain biological relevance; their statistical significance likely reflects low-level noise.
- **Key take-away:** Despite reaching nominal FDR significance, none of the eight signatures shows a large, clean separation between the two groups. Medians are close, IQRs broadly overlap, and many samples in each group are indistinguishable. The effect sizes are biologically modest.

---

### Plot 5 — Volcano Plot

![Volcano Plot: Signature Differences (Early vs Late Onset)](results/plots/plot5_volcano.jpg)

**What is shown:** A volcano plot where each point represents one COSMIC SBS signature. The x-axis shows the Log2 fold change (Early/Late), with positive values indicating higher exposure in early-onset. The y-axis shows the −log10 of the FDR-adjusted p-value, so higher positions indicate greater statistical significance. The horizontal dashed line marks the FDR = 0.05 threshold; the two vertical dashed lines flank |Log2FC| = 0.5.

**Detailed observations from this plot:**

- **No signature crosses the FDR = 0.05 significance threshold (horizontal dashed line).** The horizontal dashed line sits at the very top of the y-axis (at approximately −log10(0.05) ≈ 1.3), and all points cluster well below it. This is the central statistical finding of the analysis: after correcting for multiple comparisons across all tested COSMIC signatures, **no individual COSMIC signature is significantly differentially active between early-onset and late-onset breast cancer** at the FDR 5% level in this cohort.
- **One outlier point lies in the upper-left quadrant** (approximately Log2FC ≈ −8 to −10, −log10(FDR) ≈ 0.65). This signature has a large negative fold change — meaning it is nominally higher in late-onset — but its adjusted p-value is still far from significance. The extreme negative Log2FC is likely driven by a signature with very low expression in early-onset samples (near-zero denominator), making the fold change numerically unstable rather than biologically meaningful.
- **One outlier point lies far to the right** (approximately Log2FC ≈ 13–14, −log10(FDR) ≈ 0.13), indicating a signature nominally much higher in early-onset with an extreme fold change. Again, this almost certainly reflects a near-zero mean in the late-onset group combined with very sparse non-zero exposures in a few early-onset outliers — a statistical artefact rather than a genuine biological enrichment.
- **The majority of signatures cluster tightly near the origin** (Log2FC between −1 and +1, −log10(FDR) near 0), indicating that most COSMIC signatures have very similar exposure distributions between the two groups and contribute negligible statistical signal.
- **The absence of red labelled points** (which would denote signatures meeting both |Log2FC| > 0.5 and FDR < 0.05 thresholds simultaneously) confirms that no signature is both statistically significant and biologically meaningful by the pre-specified criteria.
- **Biological interpretation of the null result:** This null result does not imply that early- and late-onset breast cancer are biologically identical. Rather, it suggests that: (1) the effect sizes of age-related signature differences are small relative to the within-group heterogeneity; (2) WXS data may have insufficient mutation counts to reliably decompose signatures at the individual sample level; and (3) the pooled analysis across all molecular subtypes (ER+/HER2+/TNBC) may dilute subtype-specific age effects. Subtype-stratified analyses and WGS data would provide greater power.

---

### Plot 6 — Cosine Similarity (Reconstruction Quality)

Assesses how well the COSMIC signature set reconstructs each sample's observed mutation profile. Cosine similarity ranges from 0 (no match) to 1 (perfect reconstruction). Values above 0.9 indicate that the fitted COSMIC signatures adequately explain the sample's mutational landscape. Samples with low cosine similarity (<0.8) may carry novel or poorly characterised processes not represented in COSMIC v3.3.1. A bimodal distribution within a group would suggest two subpopulations with fundamentally different mutational landscapes.

---

## Biological Interpretation

### Expected findings based on literature

| Finding | Expected Direction | Biological Explanation |
|---------|-------------------|----------------------|
| SBS3 (HRD) | ↑ Early-onset | Higher prevalence of BRCA1/2 germline mutations in younger patients |
| SBS2/SBS13 (APOBEC) | ↑ Early-onset | APOBEC-driven mutagenesis more prominent in younger, hormone-driven tumours |
| SBS5 (Clock-like) | ↑ Late-onset | Accumulated replication errors over longer lifetime |
| SBS1 (CpG deamination) | ↑ Late-onset | Age-associated methylation-driven mutations |

### What the data actually show

The results reveal **no statistically significant differences in COSMIC signature exposures between early-onset and late-onset breast cancer** after FDR correction (Plot 5). Directional trends exist — SBS2 (APOBEC) is modestly higher in early-onset, and SBS1 is nominally higher in late-onset — but these do not survive multiple-testing correction. The mutation spectrum (Plot 1) shows the two groups share the same dominant mutational processes, differing only quantitatively at specific trinucleotide contexts. The heatmap (Plot 3) confirms that onset age does not stratify tumours into distinct mutational subgroups. Taken together, the findings suggest that **age at onset, when considered independently of molecular subtype, is a weak predictor of COSMIC mutational signature composition in this cohort**.

---

## Limitations

- Analysis restricted to WXS data; whole-genome sequencing would provide higher mutation counts and more reliable signature fitting
- The Mann-Whitney U test is non-parametric but does not account for tumour purity, subtype heterogeneity (e.g. ER+/HER2+/TNBC), or technical batch effects
- COSMIC refitting assumes all mutations are explained by known signatures; novel or composite processes may be misattributed
- Age cutoffs (≤45 / ≥55) exclude intermediate-age patients to maximise contrast but reduce sample size
- Low per-sample mutation burden in WXS data limits reliable multi-signature decomposition; signatures with sparse exposures may reflect refitting noise rather than genuine biological activity

---

## References

1. Alexandrov et al. (2020). The repertoire of mutational signatures in human cancer. *Nature*, 578, 94–101.
2. COSMIC Mutational Signatures v3.3.1 — https://cancer.sanger.ac.uk/signatures/
3. Blokzijl et al. (2018). MutationalPatterns: comprehensive genome-wide analysis of mutational processes. *Genome Medicine*, 10, 33.
4. Colaprico et al. (2016). TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Research*, 44(8), e71.
5. Mayakonda et al. (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research*, 28, 1747–1756.

---

## Author

**Poulami** — [@g-Poulami](https://github.com/g-Poulami)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
