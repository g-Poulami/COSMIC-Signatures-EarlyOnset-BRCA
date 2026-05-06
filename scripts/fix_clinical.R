library(jsonlite)

url <- "https://api.gdc.cancer.gov/cases?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22project.project_id%22%2C%22value%22%3A%5B%22TCGA-BRCA%22%5D%7D%7D&fields=submitter_id,diagnoses.age_at_diagnosis&format=json&size=2000"

res  <- fromJSON(url)
hits <- res$data$hits

diag_df <- lapply(seq_len(nrow(hits)), function(i) {
    age <- tryCatch({
        val <- hits$diagnoses[[i]]$age_at_diagnosis
        if (is.null(val) || length(val) == 0) NA_real_ else as.numeric(val[1])
    }, error = function(e) NA_real_)
    data.frame(submitter_id = hits$submitter_id[i],
               age_at_index = age / 365.25)
})

clinical <- do.call(rbind, diag_df)
saveRDS(clinical, "clinical.rds")
message(sprintf("Saved %d patients to clinical.rds", nrow(clinical)))
print(head(clinical))
