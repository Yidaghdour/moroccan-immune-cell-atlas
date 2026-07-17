#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(motifmatchr)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Matrix)
})

parse_args <- function() {
  defaults <- list(
    regions = "results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz",
    jaspar = "references/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
    hocomoco = "references/motifs/HOCOMOCO_v14_H14CORE_jaspar_explicit.txt",
    output = "results/tfbs_gene_overlap/tfbs_motif_hits.tsv.gz",
    `p-cutoff` = "5e-5"
  )
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- substring(key, 3L)
    if (identical(key, "help")) {
      cat("Arguments:\n")
      cat("  --regions PATH\n")
      cat("  --jaspar PATH\n")
      cat("  --hocomoco PATH\n")
      cat("  --output PATH\n")
      cat("  --p-cutoff NUMERIC\n")
      quit(save = "no", status = 0L)
    }
    if (i == length(args)) {
      stop("Missing value for --", key)
    }
    defaults[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  list(
    regions = defaults$regions,
    jaspar = defaults$jaspar,
    hocomoco = defaults$hocomoco,
    output = defaults$output,
    p_cutoff = as.numeric(defaults$`p-cutoff`)
  )
}

ensure_parent <- function(path) {
  dir.create(dirname(normalizePath(path, mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
}

write_tsv_auto <- function(dt, path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    write.table(dt, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    fwrite(dt, path, sep = "\t")
  }
}

scan_one_db <- function(db_name, motif_path, gr, p_cutoff) {
  message("Reading ", db_name, " motifs from ", motif_path)
  motifs <- readJASPARMatrix(motif_path, matrixClass = "PFM")
  motif_ids <- vapply(motifs, ID, character(1))
  motif_names <- vapply(motifs, name, character(1))

  message("Scanning ", db_name, " motifs (", length(motifs), " matrices)...")
  motif_ix <- matchMotifs(
    motifs,
    gr,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    p.cutoff = p_cutoff
  )
  motif_mat <- motifMatches(motif_ix)
  match_idx <- which(motif_mat, arr.ind = TRUE)
  if (nrow(match_idx) == 0) {
    return(data.table(
      region_id = character(),
      motif_db = character(),
      motif_id = character(),
      motif_name = character()
    ))
  }
  data.table(
    region_id = rownames(motif_mat)[match_idx[, "row"]],
    motif_db = db_name,
    motif_id = motif_ids[match_idx[, "col"]],
    motif_name = motif_names[match_idx[, "col"]]
  )
}

opt <- parse_args()

message("Loading target regions...")
regions <- fread(opt$regions)
required <- c("region_id", "chrom", "start", "end")
missing <- setdiff(required, names(regions))
if (length(missing) > 0) {
  stop("Region table is missing required columns: ", paste(missing, collapse = ", "))
}
regions <- unique(regions[, .(region_id, chrom, start, end)])
regions <- regions[end > start]
if (nrow(regions) == 0) {
  stop("No valid regions to scan.")
}

gr <- GRanges(
  seqnames = regions$chrom,
  ranges = IRanges(start = as.integer(regions$start) + 1L, end = as.integer(regions$end))
)
names(gr) <- regions$region_id

hits <- rbindlist(
  list(
    scan_one_db("JASPAR2024_CORE_vertebrates", opt$jaspar, gr, opt$p_cutoff),
    scan_one_db("HOCOMOCO_v14_H14CORE", opt$hocomoco, gr, opt$p_cutoff)
  ),
  use.names = TRUE,
  fill = TRUE
)

ensure_parent(opt$output)
write_tsv_auto(hits, opt$output)
message("Wrote ", nrow(hits), " motif hits to ", normalizePath(opt$output, mustWork = FALSE))
