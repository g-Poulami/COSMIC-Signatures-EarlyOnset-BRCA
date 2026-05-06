# 🧬 COSMIC Mutational Signatures in Early-Onset vs Late-Onset Breast Cancer

[![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.17+-1a9641?style=flat-square)](https://bioconductor.org/)
[![TCGA](https://img.shields.io/badge/Data-TCGA--BRCA-blue?style=flat-square)](https://portal.gdc.cancer.gov/)
[![COSMIC](https://img.shields.io/badge/Signatures-COSMIC_v3.3.1-orange?style=flat-square)](https://cancer.sanger.ac.uk/signatures/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)]()
[![GDC Client](https://img.shields.io/badge/GDC_Client-v1.6.1-lightgrey?style=flat-square)](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20WSL2-informational?style=flat-square&logo=linux&logoColor=white)]()

---

## 📌 Overview

This project investigates differences in **somatic mutational signatures** between **early-onset** (≤45 years) and **late-onset** (≥55 years) breast cancer patients using whole-exome sequencing (WXS) data from **The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA)** cohort.

Mutational signatures reflect the underlying biological processes — DNA repair defects, APOBEC activity, ageing-related damage — that drive tumour evolution. By comparing signature exposures across age groups, this analysis aims to uncover **age-specific mutagenic mechanisms** in breast cancer.

---

## 🧪 Scientific Background

| Signature | Biological Process | Relevance to BRCA |
|-----------|-------------------|-------------------|
| **SBS3**  | Homologous recombination deficiency (HRD) | BRCA1/2 dysfunction |
| **SBS2**  | APOBEC cytidine deaminase activity | Elevated in younger patients |
| **SBS13** | APOBEC cytidine deaminase activity | Co-occurs with SBS2 |
| **SBS5**  | Clock-like, age-associated damage | Higher in older patients |
| **SBS1**  | Spontaneous deamination of 5-methylcytosine | Age-correlated |

Early-onset breast cancer is often associated with germline BRCA1/2 mutations and HRD (SBS3), while late-onset disease tends to show more clock-like signatures (SBS5, SBS1) due to accumulated replication errors over time.

---

## 📂 Repository Structure

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
│   │   ├── plot1_mutation_spectrum.pdf    # 96-channel trinucleotide profiles
│   │   ├── plot2_stacked_bar.pdf         # Per-sample stacked bar chart
│   │   ├── plot3_heatmap.pdf             # Signature exposure heatmap
│   │   ├── plot4_boxplots.pdf            # Signature boxplots by group
│   │   ├── plot5_volcano.pdf             # Volcano plot
│   │   └── plot6_cosine_similarity.pdf   # Reconstruction quality
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

## 🔬 Data

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

## ⚙️ Dependencies

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

## 🚀 Reproducing the Analysis

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

## 📊 Outputs & Interpretation

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
**`results/plots/plot1_mutation_spectrum.pdf`**

Shows the relative frequency of all 96 single base substitution (SBS) mutation types (6 substitution classes × 16 trinucleotide contexts) averaged across early-onset and late-onset patients separately.

**How to read it:** Each bar represents one of the 96 trinucleotide mutation contexts. The 6 colours correspond to the 6 base substitution types (C>A, C>G, C>T, T>A, T>C, T>G). Differences in bar heights between groups indicate shifts in the dominant mutational processes.

**What to look for:**
- Elevated C>T transitions at TpCpN contexts → APOBEC activity (SBS2/SBS13), potentially higher in early-onset
- Elevated C>T at CpG dinucleotides → clock-like deamination (SBS1), potentially higher in late-onset

---

### Plot 2 — Stacked Bar Chart of Signature Contributions
**`results/plots/plot2_stacked_bar.pdf`**

Displays the **relative contribution** of each COSMIC signature to every individual sample, with samples split by onset group. Only signatures with a mean relative exposure >1% across all samples are shown.

**How to read it:** Each vertical bar represents one patient sample. The coloured segments show the proportion of mutations attributed to each signature. Samples are ordered within each group.

**What to look for:**
- Heterogeneity within groups (some samples dominated by SBS3, others by SBS5)
- Whether early-onset samples show more SBS3/SBS2 enrichment as a group
- Outlier samples with unusual signature profiles

---

### Plot 3 — Signature Exposure Heatmap
**`results/plots/plot3_heatmap.pdf`**

A clustered heatmap showing relative signature exposures across all samples (columns) and signatures (rows), with samples annotated by onset group (red = early, blue = late).

**How to read it:** Darker colours indicate higher relative exposure. Hierarchical clustering groups samples with similar signature profiles together. The colour annotation bar at the top shows group membership.

**What to look for:**
- Clusters of early-onset samples co-localising with high SBS3 exposure
- Separation of groups along the dendrogram suggesting distinct mutational landscapes
- Signatures that are uniformly high or low across all samples

---

### Plot 4 — Signature Boxplots by Group
**`results/plots/plot4_boxplots.pdf`**

Boxplots comparing the **relative exposure** of statistically significant signatures (FDR < 0.05) between early-onset and late-onset groups, with individual sample points overlaid as a jitter plot.

**How to read it:** Each panel shows one signature. The box spans the interquartile range (IQR), the line is the median, and whiskers extend to 1.5× IQR. Individual dots represent patients.

**What to look for:**
- Signatures where the median and IQR are clearly separated between groups
- High within-group variance (wide boxes) indicating heterogeneity
- Outlier samples with extremely high exposures

---

### Plot 5 — Volcano Plot
**`results/plots/plot5_volcano.pdf`**

A volcano plot displaying **Log2 fold change** (x-axis, Early/Late) against **−log10 FDR-adjusted p-value** (y-axis) for all tested COSMIC signatures.

**How to read it:**
- Points in the **upper right** → significantly higher in early-onset
- Points in the **upper left** → significantly higher in late-onset
- Dashed horizontal line → FDR = 0.05 significance threshold
- Dashed vertical lines → |Log2FC| = 0.5 effect size threshold
- Red labelled points → statistically and biologically significant signatures

**What to look for:**
- SBS3 in the upper right (HRD enrichment in early-onset)
- SBS5/SBS1 in the upper left (clock-like enrichment in late-onset)
- Signatures near the origin → no meaningful difference between groups

---

### Plot 6 — Cosine Similarity (Reconstruction Quality)
**`results/plots/plot6_cosine_similarity.pdf`**

A violin + boxplot showing the **cosine similarity** between each sample's observed mutational profile and the profile reconstructed from the fitted COSMIC signatures.

**How to read it:** Cosine similarity ranges from 0 (completely different) to 1 (perfect reconstruction). Higher values mean the COSMIC signatures explain the sample's mutations well. Values >0.9 are generally considered good fits.

**What to look for:**
- Samples with low cosine similarity (<0.8) may harbour novel or poorly characterised mutational processes
- Differences between groups in reconstruction quality could indicate that the COSMIC reference explains one group better than the other
- A bimodal distribution within a group suggests two subpopulations with different mutational landscapes

---

## 🧠 Biological Interpretation

### Expected findings based on literature

| Finding | Expected Direction | Biological Explanation |
|---------|-------------------|----------------------|
| SBS3 (HRD) | ↑ Early-onset | Higher prevalence of BRCA1/2 germline mutations in younger patients |
| SBS2/SBS13 (APOBEC) | ↑ Early-onset | APOBEC-driven mutagenesis more prominent in younger, hormone-driven tumours |
| SBS5 (Clock-like) | ↑ Late-onset | Accumulated replication errors over longer lifetime |
| SBS1 (CpG deamination) | ↑ Late-onset | Age-associated methylation-driven mutations |

---

## ⚠️ Limitations

- Analysis restricted to WXS data; whole-genome sequencing would provide higher mutation counts and more reliable signature fitting
- The Mann-Whitney U test is non-parametric but does not account for tumour purity, subtype heterogeneity (e.g. ER+/HER2+/TNBC), or technical batch effects
- COSMIC refitting assumes all mutations are explained by known signatures; novel or composite processes may be misattributed
- Age cutoffs (≤45 / ≥55) exclude intermediate-age patients to maximise contrast but reduce sample size

---

## 📚 References

1. Alexandrov et al. (2020). The repertoire of mutational signatures in human cancer. *Nature*, 578, 94–101.
2. COSMIC Mutational Signatures v3.3.1 — https://cancer.sanger.ac.uk/signatures/
3. Blokzijl et al. (2018). MutationalPatterns: comprehensive genome-wide analysis of mutational processes. *Genome Medicine*, 10, 33.
4. Colaprico et al. (2016). TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Research*, 44(8), e71.
5. Mayakonda et al. (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research*, 28, 1747–1756.

---

## 👩‍💻 Author

**Poulami** — [@g-Poulami](https://github.com/g-Poulami)

---

## 📄 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
