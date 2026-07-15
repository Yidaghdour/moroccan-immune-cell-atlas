library(edgeR)
library(tidyverse)

setwd("/.../.../")


# ============================================================
# FILES
# ============================================================

counts_file <- "PBMC_7major_ATAC_counts_LPS_Female.tsv"
meta_file   <- "PBMC_7major_ATAC_metadata_LPS_Female.tsv"
annot_file  <- "PBMC_7major_ATAC_peak_annotations.tsv"

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

meta$Lifestyle <- factor(
  toupper(meta$Lifestyle),
  levels = c("RURAL", "URBAN")
)

meta$Age <- as.numeric(meta$Age)

cat("Donors per group:\n")
print(table(meta$Lifestyle))

# ============================================================
# CREATE edgeR OBJECT
# ============================================================

y <- DGEList(
  counts = counts,
  samples = meta
)

design <- model.matrix(~ Lifestyle + Age, data = meta)

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

pdf("LPS_Female_edgeR_BCV_minCount5.pdf", width = 6, height = 6)
plotBCV(y)
dev.off()

fit <- glmQLFit(
  y,
  design,
  robust = TRUE
)

pdf("LPS_Female_edgeR_QLDisp_minCount5.pdf", width = 6, height = 6)
plotQLDisp(fit)
dev.off()

# ============================================================
# TEST URBAN VS RURAL
# ============================================================
# logFC > 0 = higher accessibility in Urban
# logFC < 0 = higher accessibility in Rural

qlf <- glmQLFTest(
  fit,
  coef = "LifestyleURBAN"
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
      FDR < 0.05 & logFC > 0 ~ "Urban_higher",
      FDR < 0.05 & logFC < 0 ~ "Rural_higher",
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
  "LPS_Female_edgeR_results_minCount5_WITH_BEST_GENE.csv",
  row.names = FALSE
)

sig_res_annotated <- res_annotated %>%
  filter(FDR < 0.05)

write.csv(
  sig_res_annotated,
  "LPS_Female_edgeR_significant_FDR05_minCount5_WITH_BEST_GENE.csv",
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
  "LPS_Female_edgeR_DGEList_minCount5.rds"
)

saveRDS(
  fit,
  "LPS_Female_edgeR_fit_minCount5.rds"
)

saveRDS(
  qlf,
  "LPS_Female_edgeR_qlf_minCount5.rds"
)



peak_id <- "chr17-36104295-36107029"

plot_df <- res_annotated(
  Sample = colnames(logCPM),
  logCPM = as.numeric(logCPM[peak_id, ]),
  Lifestyle = factor(
    toupper(meta$Lifestyle),
    levels = c("RURAL", "URBAN")
  )
)

# --------------------------------------------------
# Plot
# --------------------------------------------------

ggplot(res_annotated,
       aes(x = Lifestyle,
           y = logCPM,
           fill = Lifestyle)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "CCL4-associated peak accessibility",
    subtitle = peak_id,
    x = "",
    y = "logCPM (TMM normalized)"
  )












library(edgeR)
library(ggplot2)

# --------------------------------------------------
# Load data
# --------------------------------------------------

counts <- read.delim(
  "PBMC_7major_ATAC_counts_LPS_Female.tsv",
  row.names = 1,
  check.names = FALSE
)

meta <- read.delim(
  "PBMC_7major_ATAC_metadata_LPS_Female.tsv",
  check.names = FALSE
)

# Match order
meta <- meta[match(colnames(counts), meta$pb_sample_id), ]

# --------------------------------------------------
# edgeR normalization
# --------------------------------------------------

y <- DGEList(counts = counts)

keep <- filterByExpr(
  y,
  group = factor(meta$Lifestyle),
  min.count = 5
)

y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)

logCPM <- cpm(
  y,
  log = TRUE,
  prior.count = 1
)

# --------------------------------------------------
# Extract CCL4 peak
# --------------------------------------------------

peak_id <- "chr16-56624972-56626243"

plot_df <- data.frame(
  Sample = colnames(logCPM),
  logCPM = as.numeric(logCPM[peak_id, ]),
  Lifestyle = factor(
    toupper(meta$Lifestyle),
    levels = c("RURAL", "URBAN")
  )
)

# --------------------------------------------------
# Plot
# --------------------------------------------------

ggplot(plot_df,
       aes(x = Lifestyle,
           y = logCPM,
           fill = Lifestyle)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.08, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "MT1E-associated peak accessibility",
    subtitle = peak_id,
    x = "",
    y = "logCPM (TMM normalized)"
  )
