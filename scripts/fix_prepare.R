# Read all MAF files directly, handling type mismatches
library(data.table)

maf_files <- list.files(
    "GDCdata/TCGA-BRCA/Simple_Nucleotide_Variation/Masked_Somatic_Mutation",
    pattern = "\\.maf\\.gz$",
    recursive = TRUE,
    full.names = TRUE
)

message(sprintf("Reading %d MAF files...", length(maf_files)))

maf_list <- lapply(maf_files, function(f) {
    tryCatch({
        dt <- fread(f, skip = "Hugo_Symbol", showProgress = FALSE)
        # Force Tumor_Seq_Allele2 to character to avoid type conflicts
        dt[, Tumor_Seq_Allele2 := as.character(Tumor_Seq_Allele2)]
        dt
    }, error = function(e) {
        message(sprintf("Skipping %s: %s", basename(f), e$message))
        NULL
    })
})

# Remove NULLs (failed files)
maf_list <- Filter(Negate(is.null), maf_list)
message(sprintf("Successfully read %d files", length(maf_list)))

# Combine all, forcing all columns to consistent types
maf_data <- rbindlist(maf_list, fill = TRUE, use.names = TRUE)
maf_data[, Tumor_Seq_Allele2 := as.character(Tumor_Seq_Allele2)]

saveRDS(maf_data, "maf_combined.rds")
message("Saved to maf_combined.rds")
