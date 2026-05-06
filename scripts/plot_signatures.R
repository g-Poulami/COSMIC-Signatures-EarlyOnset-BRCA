# ============================================================
# Comprehensive Visualization Script
# TCGA-BRCA Mutational Signature Analysis
# Early-Onset vs Late-Onset
# ============================================================

library(MutationalPatterns)
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(TCGAbiolinks)

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading data...")
maf_data <- readRDS("maf_combined.rds")

clinical <- readRDS("clinical.rds")
early_ids <- clinical$submitter_id[!is.na(clinical$age_at_index) & clinical$age_at_index <= 45]
late_ids  <- clinical$submitter_id[!is.na(clinical$age_at_index) & clinical$age_at_index >= 55]

# ── Build mutation matrix ──────────────────────────────────────────────────────
message("Building mutation matrix...")
maf <- read.maf(maf = maf_data)

gr <- GRanges(
    seqnames = maf@data$Chromosome,
    ranges   = IRanges(start = maf@data$Start_Position,
                       end   = maf@data$End_Position),
    REF      = maf@data$Reference_Allele,
    ALT      = maf@data$Tumor_Seq_Allele2,
    sample   = maf@data$Tumor_Sample_Barcode
)
gr <- gr[nchar(gr$REF) == 1 & nchar(gr$ALT) == 1]
GenomeInfoDb::genome(gr) <- "hg38"
GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"

sample_names <- unique(gr$sample)
grl <- GRangesList(lapply(sample_names, function(s) gr[gr$sample == s]))
names(grl) <- sample_names
GenomeInfoDb::seqlevelsStyle(grl) <- "UCSC"
GenomeInfoDb::genome(grl) <- "hg38"

ref_genome <- BSgenome.Hsapiens.UCSC.hg38
mut_matrix <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

# ── COSMIC signatures ─────────────────────────────────────────────────────────
message("Loading COSMIC v3.3.1 signatures...")
cosmic_sigs <- read.table(
    "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh38.txt",
    sep = "\t", header = TRUE, row.names = 1
)
cosmic_sigs <- cosmic_sigs[rownames(mut_matrix), , drop = FALSE]
cosmic_matrix <- as.matrix(cosmic_sigs)

fit_res <- fit_to_signatures(mut_matrix, cosmic_matrix)

# ── Group labels ──────────────────────────────────────────────────────────────
exposures     <- as.data.frame(t(fit_res$contribution))
exposures$pid <- substr(rownames(exposures), 1, 12)
exposures$group <- case_when(
    exposures$pid %in% early_ids ~ "Early-onset",
    exposures$pid %in% late_ids  ~ "Late-onset",
    TRUE ~ NA_character_
)
exposures_grp <- exposures[!is.na(exposures$group), ]

# Relative exposures
rel_exp        <- sweep(fit_res$contribution, 2, colSums(fit_res$contribution), "/")
rel_df         <- as.data.frame(t(rel_exp))
rel_df$pid     <- substr(rownames(rel_df), 1, 12)
rel_df$group   <- case_when(
    rel_df$pid %in% early_ids ~ "Early-onset",
    rel_df$pid %in% late_ids  ~ "Late-onset",
    TRUE ~ NA_character_
)
rel_grp <- rel_df[!is.na(rel_df$group), ]

sig_cols <- setdiff(colnames(exposures_grp), c("pid", "group"))

# ── Stats for all signatures ───────────────────────────────────────────────────
results <- lapply(sig_cols, function(sig) {
    e <- exposures_grp[[sig]][exposures_grp$group == "Early-onset"]
    l <- exposures_grp[[sig]][exposures_grp$group == "Late-onset"]
    if (length(e) < 3 || length(l) < 3) return(NULL)
    test <- wilcox.test(e, l, exact = FALSE)
    fc   <- (median(e, na.rm=TRUE) + 1e-6) / (median(l, na.rm=TRUE) + 1e-6)
    data.frame(
        Signature    = sig,
        Early_median = median(e, na.rm=TRUE),
        Late_median  = median(l, na.rm=TRUE),
        Log2FC       = log2(fc),
        P_value      = test$p.value
    )
})
results_df       <- do.call(rbind, Filter(Negate(is.null), results))
results_df$P_adj <- p.adjust(results_df$P_value, method = "BH")
write.csv(results_df, "signature_stats_results.csv", row.names = FALSE)

GROUP_COLORS <- c("Early-onset" = "#E63946", "Late-onset" = "#457B9D")

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 1: 96-channel mutation spectrum per group
# ════════════════════════════════════════════════════════════════════════════════
message("Plot 1: Mutation spectrum...")

early_samples <- rownames(exposures)[substr(rownames(exposures), 1, 12) %in% early_ids]
late_samples  <- rownames(exposures)[substr(rownames(exposures), 1, 12) %in% late_ids]

early_profile <- rowSums(mut_matrix[, colnames(mut_matrix) %in% early_samples, drop=FALSE])
late_profile  <- rowSums(mut_matrix[, colnames(mut_matrix) %in% late_samples,  drop=FALSE])

profile_mat <- cbind(
    `Early-onset` = early_profile / sum(early_profile),
    `Late-onset`  = late_profile  / sum(late_profile)
)

pdf("plot1_mutation_spectrum.pdf", width = 14, height = 6)
plot_96_profile(profile_mat, condensed = FALSE,
                ymax = max(profile_mat) * 1.15) +
    ggtitle("96-Channel Trinucleotide Mutation Spectrum by Onset Group") +
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
dev.off()
message("  Saved plot1_mutation_spectrum.pdf")

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 2: Stacked bar chart — relative signature contributions per sample
# ════════════════════════════════════════════════════════════════════════════════
message("Plot 2: Stacked bar chart...")

# Keep only signatures with mean relative exposure > 1%
keep_sigs <- names(which(rowMeans(rel_exp) > 0.01))

stack_df <- rel_grp %>%
    select(pid, group, all_of(keep_sigs)) %>%
    pivot_longer(-c(pid, group), names_to = "Signature", values_to = "Exposure") %>%
    arrange(group, pid)

stack_df$pid <- factor(stack_df$pid,
    levels = unique(stack_df$pid[order(stack_df$group)]))

n_sigs  <- length(keep_sigs)
pal     <- colorRampPalette(brewer.pal(12, "Set3"))(n_sigs)

pdf("plot2_stacked_bar.pdf", width = 18, height = 6)
ggplot(stack_df, aes(x = pid, y = Exposure, fill = Signature)) +
    geom_bar(stat = "identity", width = 1) +
    facet_wrap(~group, scales = "free_x") +
    scale_fill_manual(values = setNames(pal, keep_sigs)) +
    labs(title = "Relative COSMIC Signature Contributions per Sample",
         x = "Sample", y = "Relative Exposure") +
    theme_bw(base_size = 10) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title   = element_text(hjust = 0.5, face = "bold"),
          legend.key.size = unit(0.4, "cm"))
dev.off()
message("  Saved plot2_stacked_bar.pdf")

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 3: Signature exposure heatmap
# ════════════════════════════════════════════════════════════════════════════════
message("Plot 3: Heatmap...")

heat_mat <- t(rel_exp[keep_sigs,
    c(early_samples[early_samples %in% colnames(rel_exp)],
      late_samples[ late_samples  %in% colnames(rel_exp)])])

ann_df <- data.frame(
    Group = ifelse(rownames(heat_mat) %in% early_samples, "Early-onset", "Late-onset"),
    row.names = rownames(heat_mat)
)
ann_colors <- list(Group = GROUP_COLORS)

pdf("plot3_heatmap.pdf", width = 12, height = 10)
pheatmap(
    t(heat_mat),
    annotation_col  = ann_df,
    annotation_colors = ann_colors,
    show_colnames   = FALSE,
    cluster_cols    = TRUE,
    cluster_rows    = TRUE,
    color           = colorRampPalette(c("white", "#FEE08B", "#D73027"))(100),
    main            = "Signature Exposure Heatmap (Early vs Late Onset)",
    fontsize        = 9,
    border_color    = NA
)
dev.off()
message("  Saved plot3_heatmap.pdf")

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 4: Boxplots for all significant signatures (FDR < 0.05)
# ════════════════════════════════════════════════════════════════════════════════
message("Plot 4: Signature boxplots...")

sig_sigs <- intersect(results_df$Signature[results_df$P_adj < 0.05], colnames(rel_grp))
if (length(sig_sigs) == 0) {
    message("  No FDR-significant signatures; showing top 8 by p-value instead")
    sig_sigs <- intersect(head(results_df$Signature[order(results_df$P_value)], 8), colnames(rel_grp))
}

box_df <- rel_grp %>%
    select(pid, group, all_of(sig_sigs)) %>%
    pivot_longer(-c(pid, group), names_to = "Signature", values_to = "Exposure")

# Add significance labels
sig_labels <- results_df %>%
    filter(Signature %in% sig_sigs) %>%
    mutate(label = case_when(
        P_adj < 0.001 ~ "***",
        P_adj < 0.01  ~ "**",
        P_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
    ))

pdf("plot4_boxplots.pdf", width = 14, height = 8)
ggplot(box_df, aes(x = group, y = Exposure, fill = group)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
    facet_wrap(~Signature, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = GROUP_COLORS) +
    labs(title = "Signature Exposures: Early vs Late Onset BRCA",
         x = NULL, y = "Relative Exposure", fill = "Group") +
    theme_bw(base_size = 10) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          strip.background = element_rect(fill = "#F0F0F0"),
          legend.position  = "bottom")
dev.off()
message("  Saved plot4_boxplots.pdf")

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 5: Volcano plot — effect size vs significance
# ════════════════════════════════════════════════════════════════════════════════
message("Plot 5: Volcano plot...")

volcano_df <- results_df %>%
    mutate(
        neg_log10_p = -log10(P_adj + 1e-10),
        Significant = P_adj < 0.05 & abs(Log2FC) > 0.5,
        label       = ifelse(Significant, Signature, "")
    )

pdf("plot5_volcano.pdf", width = 8, height = 6)
ggplot(volcano_df, aes(x = Log2FC, y = neg_log10_p,
                        color = Significant, label = label)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_text(vjust = -0.6, size = 3, fontface = "bold") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-0.5, 0.5),  linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E63946")) +
    labs(
        title    = "Volcano Plot: Signature Differences (Early vs Late Onset)",
        subtitle = "Positive Log2FC = higher in Early-onset",
        x        = "Log2 Fold Change (Early / Late)",
        y        = "-log10(FDR-adjusted p-value)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          legend.position = "none")
dev.off()
message("  Saved plot5_volcano.pdf")

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 6: Cosine similarity — reconstructed vs original
# ════════════════════════════════════════════════════════════════════════════════
message("Plot 6: Cosine similarity...")

cos_sim_df <- as.data.frame(cos_sim_matrix(mut_matrix, fit_res$reconstructed))
cos_diag   <- diag(as.matrix(cos_sim_df))  # per-sample similarity

cos_plot_df <- data.frame(
    sample     = names(cos_diag),
    cos_sim    = cos_diag,
    pid        = substr(names(cos_diag), 1, 12)
) %>%
    mutate(group = case_when(
        pid %in% early_ids ~ "Early-onset",
        pid %in% late_ids  ~ "Late-onset",
        TRUE ~ "Other"
    )) %>%
    filter(group != "Other")

pdf("plot6_cosine_similarity.pdf", width = 7, height = 5)
ggplot(cos_plot_df, aes(x = group, y = cos_sim, fill = group)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = GROUP_COLORS) +
    labs(
        title = "Cosine Similarity: Original vs Reconstructed Mutational Profile",
        x     = NULL,
        y     = "Cosine Similarity"
    ) +
    ylim(0, 1) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
dev.off()
message("  Saved plot6_cosine_similarity.pdf")

# ════════════════════════════════════════════════════════════════════════════════
# Summary
# ════════════════════════════════════════════════════════════════════════════════
message("\n✓ All plots saved:")
message("  plot1_mutation_spectrum.pdf")
message("  plot2_stacked_bar.pdf")
message("  plot3_heatmap.pdf")
message("  plot4_boxplots.pdf")
message("  plot5_volcano.pdf")
message("  plot6_cosine_similarity.pdf")
message("  signature_stats_results.csv")
