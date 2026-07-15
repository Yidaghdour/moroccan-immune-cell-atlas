library(edgeR)
library(tidyverse)

setwd("/.../.../")


# ============================================================
# FILES
# ============================================================

counts_file <- "Cytotoxic_T-NKcells_ATAC_counts_LPSvsUnstim_RURAL_Female.tsv"
meta_file   <- "Cytotoxic_T-NKcells_ATAC_metadata_LPSvsUnstim_RURAL_Female.tsv"
annot_file  <- "Cytotoxic_T-NKcells_ATAC_peak_annotations.tsv"

# ============================================================
# LOAD DATA
# ============================================================

counts <- read.delim(
  counts_file,
  row.names = 1,
  check.names = FALSE
)

meta <- read.delim(
  meta_file,
  check.names = FALSE
)

peak_annot <- read.delim(
  annot_file,
  check.names = FALSE
)

# ============================================================
# MATCH METADATA TO COUNTS
# ============================================================

meta <- meta %>%
  filter(pb_sample_id %in% colnames(counts)) %>%
  arrange(match(pb_sample_id, colnames(counts)))

counts <- counts[, meta$pb_sample_id]

stopifnot(all(colnames(counts) == meta$pb_sample_id))

# ============================================================
# PREPARE VARIABLES
# ============================================================
meta$Stimulation <- factor(
  dplyr::case_when(
    toupper(meta$Stimulation) == "UNSTIM" ~ "Unstim",
    toupper(meta$Stimulation) == "LPS" ~ "LPS",
    TRUE ~ NA_character_
  ),
  levels = c("Unstim", "LPS")
)

meta$Age <- as.numeric(meta$Age)

cat("Donors per group:\n")
print(table(meta$Stimulation))

# ============================================================
# CREATE edgeR OBJECT
# ============================================================

y <- DGEList(
  counts = counts,
  samples = meta
)

design <- model.matrix(~ Stimulation + Age, data = meta)

cat("\nDesign matrix:\n")
print(design)

# ============================================================
# FILTER PEAKS: min.count = 5
# ============================================================

keep <- filterByExpr(
  y,
  design = design,
  min.count = 5
)

cat("\nPeaks before filtering:", nrow(y), "\n")
cat("Peaks after filtering:", sum(keep), "\n")
cat("Percent kept:", round(100 * mean(keep), 2), "%\n")

y <- y[keep, , keep.lib.sizes = FALSE]

# ============================================================
# NORMALIZE
# ============================================================

y <- calcNormFactors(y, method = "TMM")

# ============================================================
# DISPERSION + MODEL FIT
# ============================================================

y <- estimateDisp(
  y,
  design,
  robust = TRUE
)

pdf("RURAL_Female_edgeR_BCV_minCount5.pdf", width = 6, height = 6)
plotBCV(y)
dev.off()

fit <- glmQLFit(
  y,
  design,
  robust = TRUE
)

pdf("RURAL_Female_edgeR_QLDisp_minCount5.pdf", width = 6, height = 6)
plotQLDisp(fit)
dev.off()

# ============================================================
# TEST LPS VS Unstim
# ============================================================
# logFC > 0 = higher accessibility in LPS
# logFC < 0 = higher accessibility in Unstim

qlf <- glmQLFTest(
  fit,
  coef = "StimulationLPS"
)

res <- topTags(
  qlf,
  n = Inf,
  sort.by = "PValue"
)$table

res <- res %>%
  rownames_to_column("peak_id") %>%
  mutate(
    FDR = p.adjust(PValue, method = "BH"),
    direction = case_when(
      FDR < 0.05 & logFC > 0 ~ "LPS_higher",
      FDR < 0.05 & logFC < 0 ~ "Unstim_higher",
      TRUE ~ "Not_significant"
    )
  )

# ============================================================
# ADD PEAK ANNOTATIONS / BEST GENE
# ============================================================

if (!"peak_id" %in% colnames(peak_annot)) {
  stop("peak_annot must contain a peak_id column.")
}

if (!"best_gene" %in% colnames(peak_annot)) {
  warning("No best_gene column found in annotation file.")
  peak_annot$best_gene <- NA
}

res_annotated <- res %>%
  left_join(
    peak_annot,
    by = "peak_id"
  )

# Put important columns first
priority_cols <- c(
  "peak_id",
  "best_gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR",
  "direction"
)

other_cols <- setdiff(colnames(res_annotated), priority_cols)

res_annotated <- res_annotated %>%
  select(any_of(priority_cols), all_of(other_cols))

# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(
  res_annotated,
  "RURAL_Female_edgeR_results_minCount5_WITH_BEST_GENE.csv",
  row.names = FALSE
)

sig_res_annotated <- res_annotated %>%
  filter(FDR < 0.05)

write.csv(
  sig_res_annotated,
  "RURAL_Female_edgeR_significant_FDR05_minCount5_WITH_BEST_GENE.csv",
  row.names = FALSE
)

cat("\nSignificant peaks at FDR < 0.05:\n")
print(table(sig_res_annotated$direction))

cat("\nTop significant annotated peaks:\n")
print(head(sig_res_annotated, 20))

# ============================================================
# SAVE OBJECTS
# ============================================================

saveRDS(
  y,
  "RURAL_Female_edgeR_DGEList_minCount5.rds"
)

saveRDS(
  fit,
  "RURAL_Female_edgeR_fit_minCount5.rds"
)

saveRDS(
  qlf,
  "RURAL_Female_edgeR_qlf_minCount5.rds"
)

