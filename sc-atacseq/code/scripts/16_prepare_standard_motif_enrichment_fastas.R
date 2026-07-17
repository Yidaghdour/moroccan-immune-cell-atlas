#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

parse_args <- function() {
  defaults <- list(
    target_regions = "results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz",
    promoter_background = "results/tfbs_gene_overlap_integrated/background_reference_promoters.tsv.gz",
    linked_background = "results/tfbs_gene_overlap_integrated/background_linked_peaks.tsv.gz",
    output_dir = "results/tfbs_gene_overlap_integrated/standard_motif_enrichment/fastas"
  )
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- substring(key, 3L)
    if (i == length(args)) {
      stop("Missing value for --", key)
    }
    defaults[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  defaults
}

allowed_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
region_order <- c("promoter", "tss_proximal_peak", "proximal_peak", "distal_peak")

safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) > 0, x, "NA")
}

dedupe_regions <- function(dt, prefix, region_class = NULL) {
  keep <- copy(dt)
  keep <- keep[chrom %in% allowed_chroms & end > start]
  if (!is.null(region_class)) {
    keep[, region_class := region_class]
  }
  keep <- unique(keep, by = c("chrom", "start", "end", "region_class"))
  keep[, seq_id := paste(prefix, safe_id(region_class), chrom, start, end, sep = "__")]
  keep[]
}

write_fasta <- function(dt, path) {
  if (nrow(dt) == 0L) {
    stop("No regions available for FASTA: ", path)
  }
  gr <- GRanges(
    seqnames = dt$chrom,
    ranges = IRanges(start = as.integer(dt$start) + 1L, end = as.integer(dt$end))
  )
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  names(seqs) <- dt$seq_id
  writeXStringSet(seqs, filepath = path, format = "fasta", width = 80)
}

opt <- parse_args()
out_dir <- normalizePath(opt$output_dir, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target <- fread(opt$target_regions)
promoter_bg <- fread(opt$promoter_background)
linked_bg <- fread(opt$linked_background)

target <- dedupe_regions(target, "TARGET")
promoter_bg <- dedupe_regions(promoter_bg, "BG", region_class = "promoter")
linked_bg <- dedupe_regions(linked_bg, "BG")

contrast_defs <- list(
  list(name = "promoter_all", target_classes = "promoter", background = promoter_bg),
  list(name = "linked_pooled_all", target_classes = c("tss_proximal_peak", "proximal_peak", "distal_peak"), background = linked_bg),
  list(name = "tss_proximal_peak_all", target_classes = "tss_proximal_peak", background = linked_bg[region_class == "tss_proximal_peak"]),
  list(name = "proximal_peak_all", target_classes = "proximal_peak", background = linked_bg[region_class == "proximal_peak"]),
  list(name = "distal_peak_all", target_classes = "distal_peak", background = linked_bg[region_class == "distal_peak"])
)

summary_rows <- list()
for (contrast in contrast_defs) {
  foreground <- target[region_class %in% contrast$target_classes]
  background <- contrast$background
  fg_path <- file.path(out_dir, paste0(contrast$name, ".target.fa"))
  bg_path <- file.path(out_dir, paste0(contrast$name, ".background.fa"))
  write_fasta(foreground, fg_path)
  write_fasta(background, bg_path)
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    contrast = contrast$name,
    target_region_classes = paste(contrast$target_classes, collapse = ";"),
    target_sequences = nrow(foreground),
    background_sequences = nrow(background),
    target_fasta = fg_path,
    background_fasta = bg_path
  )
}

summary <- rbindlist(summary_rows)
fwrite(summary, file.path(dirname(out_dir), "motif_enrichment_contrasts.tsv"), sep = "\t")
print(summary)
