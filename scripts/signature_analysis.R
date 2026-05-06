# Mutational Signature Analysis Pipeline
# Dataset: TCGA-BRCA
# Comparison: Early-Onset (<=45) vs Late-Onset (>=55)

library(TCGAbiolinks)
library(maftools)
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(gridExtra)

# ── 1. Data Download ──────────────────────────────────────────────────────────
message("Querying TCGA-BRCA somatic mutation data...")
query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

# FIX 1: client argument is just the path to the binary
#GDCdownload(query, method = "client")  # already downloaded
maf_data <- readRDS("maf_combined.rds")

# ── 2. Clinical Stratification ────────────────────────────────────────────────
message("Stratifying patients by age at diagnosis...")
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# Extract TCGA patient barcodes (first 12 characters of sample ID)
early_ids <- clinical$submitter_id[
    !is.na(clinical$age_at_index) & clinical$age_at_index <= 45
]
late_ids <- clinical$submitter_id[
    !is.na(clinical$age_at_index) & clinical$age_at_index >= 55
]

message(sprintf("Early-onset: %d patients | Late-onset: %d patients",
                length(early_ids), length(late_ids)))

# ── 3. Mutational Matrix (96-channel SBS) ─────────────────────────────────────
message("Building mutational matrix...")

# Load MAF
maf <- read.maf(maf = maf_data)

# FIX 2: Use MutationalPatterns throughout (consistent with fit_to_signatures)
# Convert MAF to GRanges
gr <- GRanges(
    seqnames = paste0("chr", maf@data$Chromosome),
    ranges   = IRanges(start = maf@data$Start_Position,
                       end   = maf@data$End_Position),
    REF      = maf@data$Reference_Allele,
    ALT      = maf@data$Tumor_Seq_Allele2,
    sample   = maf@data$Tumor_Sample_Barcode
)

# Keep only SNVs (MutationalPatterns requires single-base substitutions)
gr <- gr[nchar(gr$REF) == 1 & nchar(gr$ALT) == 1]

# Get reference genome
ref_genome <- BSgenome.Hsapiens.UCSC.hg38

# Build 96-channel mutation matrix
# split GRanges by sample into a GRangesList
sample_names <- unique(gr$sample)
grl <- GRangesList(lapply(sample_names, function(s) gr[gr$sample == s]))
names(grl) <- sample_names

mut_matrix <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

# ── 4. COSMIC Signature Refitting ─────────────────────────────────────────────
message("Loading COSMIC v3.3.1 signatures...")
signatures_url <- "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh38.txt"
cosmic_sigs <- read.table(signatures_url, sep = "\t", header = TRUE, row.names = 1)

# FIX 3: Fix row names — COSMIC uses e.g. "A[C>A]A", read.table converts
# brackets/> to dots; restore original format
rownames(cosmic_sigs) <- gsub("\\.", ">",
    gsub("^(.)\\.(.)\\.(.)\\.(.)$", "\\1[\\2>\\3]\\4",
         rownames(cosmic_sigs)))

# Align mutation matrix rows to COSMIC row order
cosmic_sigs <- cosmic_sigs[rownames(mut_matrix), , drop = FALSE]
cosmic_matrix <- as.matrix(cosmic_sigs)

# Refit mutations to COSMIC signatures
fit_res <- fit_to_signatures(mut_matrix, cosmic_matrix)

# ── 5. Statistical Comparison ─────────────────────────────────────────────────
message("Comparing signature exposures between groups...")

exposures <- as.data.frame(t(fit_res$contribution))  # samples x signatures
exposures$patient_id <- substr(rownames(exposures), 1, 12)

# Label groups
exposures$group <- case_when(
    exposures$patient_id %in% early_ids ~ "Early-onset",
    exposures$patient_id %in% late_ids  ~ "Late-onset",
    TRUE ~ NA_character_
)
exposures <- exposures[!is.na(exposures$group), ]

# Signatures of interest
sigs_of_interest <- c("SBS3",   # HRD / BRCA1/2
                       "SBS5",   # Clock-like
                       "SBS13",  # APOBEC
                       "SBS2")   # APOBEC

results <- lapply(sigs_of_interest, function(sig) {
    early_vals <- exposures[[sig]][exposures$group == "Early-onset"]
    late_vals  <- exposures[[sig]][exposures$group == "Late-onset"]
    test <- wilcox.test(early_vals, late_vals, exact = FALSE)
    data.frame(
        Signature   = sig,
        Early_median = median(early_vals, na.rm = TRUE),
        Late_median  = median(late_vals,  na.rm = TRUE),
        W_statistic  = test$statistic,
        P_value      = test$p.value
    )
})

results_df <- do.call(rbind, results)
results_df$P_adj <- p.adjust(results_df$P_value, method = "BH")  # FDR correction
print(results_df)

# ── 6. Visualisation ──────────────────────────────────────────────────────────
message("Plotting...")

# Relative exposures (normalise per sample)
rel_exp <- sweep(fit_res$contribution, 2, colSums(fit_res$contribution), "/")
rel_df  <- as.data.frame(t(rel_exp))
rel_df$patient_id <- substr(rownames(rel_df), 1, 12)
rel_df$group <- case_when(
    rel_df$patient_id %in% early_ids ~ "Early-onset",
    rel_df$patient_id %in% late_ids  ~ "Late-onset",
    TRUE ~ NA_character_
)
rel_df <- rel_df[!is.na(rel_df$group), ]

plot_list <- lapply(sigs_of_interest, function(sig) {
    ggplot(rel_df, aes(x = group, y = .data[[sig]], fill = group)) +
        geom_boxplot(outlier.shape = 21, alpha = 0.8) +
        geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
        scale_fill_manual(values = c("Early-onset" = "#E63946",
                                     "Late-onset"  = "#457B9D")) +
        labs(title = sig, x = NULL,
             y = "Relative exposure") +
        theme_bw(base_size = 11) +
        theme(legend.position = "none")
})

pdf("signature_comparison_BRCA.pdf", width = 10, height = 8)
grid.arrange(grobs = plot_list, ncol = 2,
             top = "COSMIC Signature Exposures: Early vs Late Onset BRCA")
dev.off()

message("Done! Results saved to signature_comparison_BRCA.pdf")
write.csv(results_df, "signature_stats_results.csv", row.names = FALSE)