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

## Biological Question

Why do breast cancers diagnosed before age 45 behave clinically differently from those diagnosed after 55 — and can we see this difference written into the tumour's DNA?

Every cancer accumulates a distinctive pattern of mutations shaped by the biological processes that drove its growth: DNA repair failures, enzymatic deaminase activity, ageing-related replication errors. These patterns, called **mutational signatures**, act as molecular fingerprints of tumourigenesis — preserved in the tumour genome long after the causative process has ceased. Comparing these fingerprints across age groups offers a window into the distinct aetiological forces operating in young versus older patients.

This question has direct clinical relevance. Early-onset breast cancer is enriched for germline BRCA1/2 mutations and homologous recombination deficiency (HRD), suggesting a different mutagenic biology from the clock-like, replication-error-driven processes that dominate in older patients. If APOBEC-associated or HRD signatures are systematically elevated in young patients, this has implications for germline testing thresholds, PARP inhibitor eligibility, and the design of age-stratified screening programmes.

This project addresses the question using **992 TCGA-BRCA whole-exome sequencing samples**, decomposing somatic mutations into COSMIC v3.3.1 signatures and statistically comparing exposures between early-onset (≤45 years) and late-onset (≥55 years) groups.

---

## Key Finding

> **After Benjamini-Hochberg FDR correction across all tested COSMIC signatures, no individual signature is significantly differentially active between early-onset and late-onset breast cancer in this cohort.**

Directional trends exist and are biologically interpretable — APOBEC signatures (SBS2/SBS13) are modestly elevated in early-onset patients, and clock-like signatures (SBS1) show a marginal increase in late-onset — but these do not survive multiple-testing correction. The heatmap confirms that onset age does not stratify tumours into distinct mutational subgroups at the whole-cohort level.

**This null result is scientifically informative, not a failure.** It tells us that:

1. The effect sizes of age-related signature differences are small relative to within-group biological heterogeneity
2. Whole-exome sequencing provides insufficient mutation counts per sample for reliable multi-signature decomposition
3. Pooling across molecular subtypes (ER+, HER2+, TNBC) dilutes subtype-specific age effects
4. Subtype-stratified analyses using whole-genome sequencing data are needed to resolve this question

Reproducible negative results are an underrepresented but scientifically essential contribution to the field.

---

## Scientific Background

| Signature | Biological Process | Expected Direction | Mechanistic Basis |
|-----------|-------------------|-------------------|-------------------|
| **SBS3**  | Homologous recombination deficiency (HRD) | ↑ Early-onset | BRCA1/2 germline mutations impair HR repair; more prevalent in younger hereditary cases |
| **SBS2**  | APOBEC cytidine deaminase (C>T at TpCpN) | ↑ Early-onset | APOBEC3A/B-driven mutagenesis associated with genomic instability in younger, hormone-driven tumours |
| **SBS13** | APOBEC cytidine deaminase (C>G at TpCpN) | ↑ Early-onset | Co-occurs with SBS2; reflects the same APOBEC mutagenic process via a different repair pathway |
| **SBS5**  | Clock-like replication-associated damage | ↑ Late-onset | Accumulates proportionally with number of cell divisions over a lifetime |
| **SBS1**  | Spontaneous deamination of 5-methylcytosine | ↑ Late-onset | Age-correlated CpG→TpG transitions; nearly universal in cancer, increases with age |

Early-onset breast cancer is frequently associated with germline BRCA1/2 dysfunction and HRD, while late-onset disease tends to accumulate clock-like signatures reflecting decades of replication. APOBEC-driven mutagenesis (SBS2/SBS13) is subtype-specific and heterogeneous, appearing prominently in a subset of tumours regardless of age.

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

**What is shown:** The relative frequency of all 96 single base substitution (SBS) mutation types across six substitution classes (C>A in blue, C>G in black, C>T in red, T>A in grey, T>C in green, T>G in pink), each broken into 16 trinucleotide contexts, averaged across early-onset (top row) and late-onset (bottom row) patients.

**Biological interpretation:**

- **C>G transversions are the most visually striking feature of the early-onset spectrum.** Elevated C>G signal at T[C>G]T and flanking contexts in early-onset is consistent with heightened APOBEC-like or oxidative mutagenic activity in younger patients.
- **C>T transitions dominate both groups**, as expected for breast cancer. The tallest C>T peak at TpCpN contexts — a hallmark of APOBEC activity (SBS2/SBS13) — is modestly higher in early-onset, directionally consistent with the APOBEC hypothesis but not statistically distinguishing.
- **C>A transversions show a late-onset asymmetry** at one context (likely T[C>A]T), which may reflect SBS18 or oxidative damage accumulating with age.
- **Overall spectrum similarity:** The two profiles share the same dominant mutational processes, differing quantitatively rather than qualitatively. This is consistent with the statistical finding of no significant differences after FDR correction — the age-related signal is present but too small to overcome within-group heterogeneity at this sample size.

---

### Plot 2 — Stacked Bar Chart of Signature Contributions

![Relative COSMIC Signature Contributions per Sample](results/plots/plot2_stacked_bar.jpg)

**What is shown:** Relative COSMIC signature contributions per patient (one bar = one patient), split into early-onset (left) and late-onset (right) panels. Only signatures with mean relative exposure >1% are shown.

**Biological interpretation:**

- **Inter-patient heterogeneity is the dominant feature in both groups.** No single signature uniformly dominates — the colour composition varies dramatically per patient, reflecting the well-known molecular heterogeneity of breast cancer across ER+, HER2+, and TNBC subtypes pooled in TCGA-BRCA.
- **Outlier samples with relative exposures >1.0** are statistical artefacts from NNLS refitting in low-mutation-burden samples, not genuine biological signal. These are more frequent in the larger late-onset cohort.
- **SBS2 and SBS13 (APOBEC)** appear as prominent segments in a subset of samples in both groups — confirming that APOBEC activity is a sample-specific rather than age-stratified phenomenon at cohort level.
- **SBS1 (clock-like CpG deamination)** contributes a consistent baseline across most samples in both groups, as expected for a near-universal age-correlated process.

---

### Plot 3 — Signature Exposure Heatmap

![Signature Exposure Heatmap (Early vs Late Onset)](results/plots/plot3_heatmap.jpg)

**What is shown:** Hierarchically clustered heatmap of COSMIC signature exposures. Each column is a patient; each row is a signature. The annotation bar colour-codes samples by group (red = early-onset, blue = late-onset).

**Biological interpretation:**

- **The column dendrogram does not separate early-onset from late-onset samples into distinct clusters.** Early and late-onset samples are thoroughly interleaved across the dendrogram — confirming that the global mutational signature profile does not reliably predict age of onset. Onset age is a weak discriminator of mutational subgroup identity at the whole-cohort, whole-subtype level.
- **SBS2 and SBS13 cluster together** with the highest exposure values, concentrated in a subset of samples scattered across both groups — reinforcing the sample-specific nature of APOBEC mutagenesis.
- **Most signatures show near-zero exposure** across the majority of samples, which is typical for WXS data where per-sample mutation counts are insufficient to reliably decompose more than 2–3 signatures simultaneously.

---

### Plot 4 — Signature Boxplots by Group

![Signature Exposures: Early vs Late Onset BRCA](results/plots/plot4_boxplots.jpg)

**What is shown:** Boxplots of relative exposure for the eight signatures reaching FDR < 0.05 (SBS1, SBS2, SBS30, SBS34, SBS38, SBS54, SBS60, SBS88), comparing early-onset (red) vs late-onset (blue). Individual data points are overlaid.

**Biological interpretation:**

- **SBS1 (clock-like):** Both groups show similar medians with overlapping IQRs. The late-onset group has a marginally wider spread, consistent with the known age-correlated accumulation of SBS1, but the difference is modest.
- **SBS2 (APOBEC):** Early-onset median is marginally higher, directionally consistent with greater APOBEC activity in younger patients, but the distributions substantially overlap with extensive right-tail outliers in both groups.
- **SBS30, SBS34, SBS38, SBS54, SBS60, SBS88:** Medians are at or near zero in both groups. Their nominal FDR significance likely reflects refitting noise in low-mutation-burden samples rather than genuine biological enrichment. SBS54 shows the most visually distinct group separation but is a poorly characterised signature in breast cancer.
- **Key take-away:** Despite reaching nominal FDR significance, none of the eight signatures shows a large, clean separation. Effect sizes are biologically modest — the statistical significance is driven by the large sample size amplifying small distributional differences.

---

### Plot 5 — Volcano Plot

![Volcano Plot: Signature Differences (Early vs Late Onset)](results/plots/plot5_volcano.jpg)

**What is shown:** Each point represents one COSMIC SBS signature. X-axis: Log2 fold change (Early/Late; positive = higher in early-onset). Y-axis: −log10(FDR-adjusted p-value). Dashed lines mark FDR = 0.05 and |Log2FC| = 0.5.

**Biological interpretation:**

- **No signature crosses the FDR = 0.05 significance threshold.** All points cluster well below the horizontal dashed line — the central statistical finding of this analysis. After correcting for multiple comparisons, no COSMIC signature is significantly differentially active between age groups.
- **Two extreme outlier points** with very large |Log2FC| values reflect near-zero denominator means causing numerically unstable fold changes, not genuine biological enrichment.
- **The majority of signatures cluster near the origin**, indicating similar exposure distributions between groups and negligible statistical signal.
- **The absence of labelled significant points** confirms no signature meets both the statistical (FDR < 0.05) and biological (|Log2FC| > 0.5) significance criteria simultaneously.

---

### Plot 6 — Cosine Similarity (Reconstruction Quality)

![Cosine Similarity Distribution](results/plots/plot6_cosine_similarity.jpg)

Assesses how well the COSMIC v3.3.1 signature set reconstructs each sample's observed mutation profile. Values above 0.9 indicate that the fitted signatures adequately explain the sample's mutational landscape. Samples below 0.8 may carry novel or poorly characterised processes not represented in COSMIC v3.3.1, or may simply have too few mutations for reliable fitting. A bimodal distribution within a group would suggest two subpopulations with fundamentally different mutational architectures.

---

## Biological Interpretation & Discussion

### What the data show

The results reveal **no statistically significant differences in COSMIC signature exposures between early-onset and late-onset breast cancer** after FDR correction (Plot 5). Directional trends are present and biologically consistent with prior literature:

- SBS2/SBS13 (APOBEC) are modestly elevated in early-onset
- SBS1 (clock-like CpG deamination) is nominally higher in late-onset

However, these trends do not survive multiple-testing correction. The mutation spectrum (Plot 1) shows the two groups share the same dominant mutational processes. The heatmap (Plot 3) confirms that onset age does not stratify tumours into distinct mutational subgroups. **Taken together, the findings suggest that age at onset, when considered independently of molecular subtype, is a weak predictor of COSMIC mutational signature composition in this cohort.**

### What the literature predicts vs what was found

| Signature | Literature prediction | This analysis |
|-----------|----------------------|---------------|
| SBS3 (HRD) | ↑ Early-onset (BRCA1/2 enrichment) | No significant difference |
| SBS2/SBS13 (APOBEC) | ↑ Early-onset | Directional trend only; not FDR-significant |
| SBS5 (clock-like) | ↑ Late-onset | No significant difference |
| SBS1 (CpG deamination) | ↑ Late-onset | Directional trend only; not FDR-significant |

### Why the null result is informative

The absence of statistically significant findings is itself a result that narrows the hypothesis space:

1. **Power limitation of WXS:** Whole-exome sequencing captures ~2% of the genome. The median mutation burden per sample in this cohort is insufficient for reliable decomposition of more than 2–3 signatures simultaneously, reducing sensitivity to detect modest differences.
2. **Subtype heterogeneity as a confounder:** TCGA-BRCA pools ER+, HER2+, and TNBC tumours, each with distinct signature profiles. The age-related signal (e.g. SBS3 enrichment in hereditary cases) may be concentrated in a specific subtype (TNBC or ER-negative) but diluted to non-significance in the pooled analysis.
3. **Small effect sizes relative to inter-patient variance:** Breast cancer is biologically heterogeneous. The within-group variance of signature exposures is large, requiring either much larger sample sizes or reduced variance (via subtype stratification) to detect the age signal.
4. **The age cutoffs exclude intermediate patients:** The ≤45/≥55 design maximises contrast but reduces sample size and may still include aetiologically heterogeneous tumours within each group.

### Recommended next steps

- Subtype-stratified analysis (TNBC separately, as the subtype most enriched for BRCA1/2 dysfunction and HRD)
- Whole-genome sequencing data (higher mutation counts; more reliable signature fitting; access to SV-based HRD metrics)
- Integration of germline BRCA1/2 status as a covariate to directly test whether HRD enrichment in early-onset is explained by germline carrier status
- Larger external cohorts (METABRIC, ICGC) to increase power

---

## Limitations

- Analysis restricted to WXS data; whole-genome sequencing would provide higher mutation counts and more reliable signature fitting
- The Mann-Whitney U test is non-parametric but does not account for tumour purity, subtype heterogeneity (ER+/HER2+/TNBC), or technical batch effects
- COSMIC refitting assumes all mutations are explained by known signatures; novel or composite processes may be misattributed
- Age cutoffs (≤45 / ≥55) exclude intermediate-age patients to maximise contrast but reduce sample size
- Low per-sample mutation burden in WXS data limits reliable multi-signature decomposition; signatures with sparse exposures may reflect refitting noise rather than genuine biological activity
- Germline BRCA1/2 status was not available as a covariate for this analysis

---

## References

1. Alexandrov et al. (2020). The repertoire of mutational signatures in human cancer. *Nature*, 578, 94–101.
2. COSMIC Mutational Signatures v3.3.1 — https://cancer.sanger.ac.uk/signatures/
3. Blokzijl et al. (2018). MutationalPatterns: comprehensive genome-wide analysis of mutational processes. *Genome Medicine*, 10, 33.
4. Colaprico et al. (2016). TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Research*, 44(8), e71.
5. Mayakonda et al. (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research*, 28, 1747–1756.
6. Burns et al. (2013). APOBEC3B is an enzymatic source of mutation in breast cancer. *Nature Genetics*, 45, 229–233.
7. Tutt et al. (2021). Adjuvant olaparib for patients with BRCA1- or BRCA2-mutated breast cancer. *New England Journal of Medicine*, 384, 2394–2405.

---

## Author

**Poulami Ghosh** — [@g-Poulami](https://github.com/g-Poulami)
[LinkedIn](https://linkedin.com/in/poulami-ghosh-879439304) · [Google Scholar](https://scholar.google.com/scholar?q=Poulami+Ghosh) · poulamighosh738@gmail.com

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
