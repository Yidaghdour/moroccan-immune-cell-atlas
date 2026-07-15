library(tidyverse)

setwd("/.../.../")

CELLTAG <- "Monocytes"

OUTDIR <- "ATAC_DID_from_existing_edgeR_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

FDR_CUTOFF <- 0.05
DID_LFC_CUTOFF <- 0.585

# ============================================================
# FUNCTION TO READ EXISTING edgeR RESULTS
# ============================================================

read_existing_atac_result <- function(lifestyle, sex) {
  
  file <- paste0(
    lifestyle, "_", sex,
    "_edgeR_results_minCount5_WITH_BEST_GENE.csv"
  )
  
  if (!file.exists(file)) {
    stop("Missing file: ", file)
  }
  
  read.csv(file, check.names = FALSE) %>%
    mutate(
      peak_id = as.character(peak_id),
      best_gene = as.character(best_gene),
      lifestyle = lifestyle,
      sex = sex
    ) %>%
    select(
      peak_id,
      best_gene,
      logFC,
      logCPM,
      F,
      PValue,
      FDR,
      direction,
      everything()
    )
}

# ============================================================
# FUNCTION TO COMPUTE DID PER SEX
# ============================================================

make_atac_did_from_existing <- function(sex_keep) {
  
  cat("\n============================================================\n")
  cat("Computing ATAC DID from existing edgeR results for:", sex_keep, "\n")
  cat("============================================================\n")
  
  urban <- read_existing_atac_result("URBAN", sex_keep) %>%
    rename(
      urban_logFC_LPS_vs_Unstim = logFC,
      urban_logCPM = logCPM,
      urban_F = F,
      urban_PValue = PValue,
      urban_FDR = FDR,
      urban_direction = direction
    ) %>%
    select(
      peak_id,
      best_gene,
      starts_with("urban_"),
      everything()
    )
  
  rural <- read_existing_atac_result("RURAL", sex_keep) %>%
    rename(
      rural_logFC_LPS_vs_Unstim = logFC,
      rural_logCPM = logCPM,
      rural_F = F,
      rural_PValue = PValue,
      rural_FDR = FDR,
      rural_direction = direction
    ) %>%
    select(
      peak_id,
      best_gene,
      starts_with("rural_"),
      everything()
    )
  
  did <- urban %>%
    inner_join(
      rural,
      by = c("peak_id", "best_gene"),
      suffix = c("_urbanfile", "_ruralfile")
    ) %>%
    mutate(
      cell_type = CELLTAG,
      sex = sex_keep,
      
      did_logFC_urban_minus_rural =
        urban_logFC_LPS_vs_Unstim - rural_logFC_LPS_vs_Unstim,
      
      response_sig_in_urban =
        urban_FDR < FDR_CUTOFF,
      
      response_sig_in_rural =
        rural_FDR < FDR_CUTOFF,
      
      response_sig_in_both =
        urban_FDR < FDR_CUTOFF &
        rural_FDR < FDR_CUTOFF,
      
      did_abs_logFC_pass =
        abs(did_logFC_urban_minus_rural) > DID_LFC_CUTOFF,
      
      did_direction = case_when(
        did_logFC_urban_minus_rural > DID_LFC_CUTOFF ~
          "Urban_stronger_LPS_opening",
        
        did_logFC_urban_minus_rural < -DID_LFC_CUTOFF ~
          "Rural_stronger_LPS_opening",
        
        TRUE ~ "No_strong_DID"
      ),
      
      heatmap_candidate =
        response_sig_in_both &
        did_abs_logFC_pass
    )
  
  # ==========================================================
  # PEAK-LEVEL RESULT
  # ==========================================================
  
  priority_cols <- c(
    "cell_type",
    "sex",
    "peak_id",
    "best_gene",
    
    "urban_logFC_LPS_vs_Unstim",
    "urban_FDR",
    "urban_direction",
    
    "rural_logFC_LPS_vs_Unstim",
    "rural_FDR",
    "rural_direction",
    
    "did_logFC_urban_minus_rural",
    "did_direction",
    
    "response_sig_in_urban",
    "response_sig_in_rural",
    "response_sig_in_both",
    "did_abs_logFC_pass",
    "heatmap_candidate"
  )
  
  did <- did %>%
    select(any_of(priority_cols), everything())
  
  write.csv(
    did,
    file.path(
      OUTDIR,
      paste0(sex_keep, "_ATAC_DID_from_existing_edgeR_peak_level.csv")
    ),
    row.names = FALSE
  )
  
  # ==========================================================
  # BEST-GENE LEVEL RESULT
  # Keep strongest DID peak per best_gene
  # ==========================================================
  
  did_gene <- did %>%
    filter(!is.na(best_gene), best_gene != "", best_gene != "nan") %>%
    group_by(best_gene) %>%
    slice_max(
      order_by = abs(did_logFC_urban_minus_rural),
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  write.csv(
    did_gene,
    file.path(
      OUTDIR,
      paste0(sex_keep, "_ATAC_DID_from_existing_edgeR_bestGene_level.csv")
    ),
    row.names = FALSE
  )
  
  did_gene_sig <- did_gene %>%
    filter(heatmap_candidate)
  
  write.csv(
    did_gene_sig,
    file.path(
      OUTDIR,
      paste0(sex_keep, "_ATAC_DID_from_existing_edgeR_heatmap_candidates_bestGene_level.csv")
    ),
    row.names = FALSE
  )
  
  cat("\nPeak-level DID direction counts:\n")
  print(table(did$did_direction))
  
  cat("\nGene-level heatmap candidate counts:\n")
  print(table(did_gene_sig$did_direction))
  
  return(did_gene)
}

# ============================================================
# RUN MALE AND FEMALE
# ============================================================

female_did <- make_atac_did_from_existing("Female")
male_did   <- make_atac_did_from_existing("Male")

combined_did <- bind_rows(female_did, male_did)

write.csv(
  combined_did,
  file.path(
    OUTDIR,
    paste0(CELLTAG, "_Male_Female_ATAC_DID_from_existing_edgeR_bestGene_level.csv")
  ),
  row.names = FALSE
)

cat("\nDone.\n")
print(dim(combined_did))

