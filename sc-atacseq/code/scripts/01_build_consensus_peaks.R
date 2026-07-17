#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Matrix)
  library(Signac)
})

default_libraries <- sprintf("1_%03d", 1:17)

parse_bool <- function(value) {
  normalized <- tolower(trimws(as.character(value)))
  if (normalized %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (normalized %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop("Expected a boolean value, received: ", value)
}

parse_args <- function() {
  defaults <- list(
    `cellranger-root` = "data/raw/atac_cellranger",
    `output-dir` = "data/processed/atac_requantified",
    libraries = paste(default_libraries, collapse = ","),
    `merge-distance` = "50",
    `bed-starts-zero-based` = "false"
  )
  args <- commandArgs(trailingOnly = TRUE)
  index <- 1L
  while (index <= length(args)) {
    key <- args[[index]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- substring(key, 3L)
    if (identical(key, "help")) {
      cat("Construct and requantify the 17-library consensus ATAC peak set.\n\n")
      cat("Arguments:\n")
      cat("  --cellranger-root PATH          Root containing LIBRARY/outs directories\n")
      cat("  --output-dir PATH               Destination for consensus and count files\n")
      cat("  --libraries ID1,ID2,...         Comma-separated library identifiers\n")
      cat("  --merge-distance INTEGER        Maximum bedtools merge gap (default: 50)\n")
      cat("  --bed-starts-zero-based BOOLEAN Treat BED starts as zero-based in GRanges\n")
      quit(save = "no", status = 0L)
    }
    if (index == length(args)) {
      stop("Missing value for --", key)
    }
    if (!key %in% names(defaults)) {
      stop("Unknown argument: --", key)
    }
    defaults[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  libraries <- trimws(strsplit(defaults$libraries, ",", fixed = TRUE)[[1]])
  libraries <- libraries[nzchar(libraries)]
  if (length(libraries) == 0L) {
    stop("At least one library must be supplied with --libraries")
  }
  list(
    cellranger_root = defaults$`cellranger-root`,
    output_dir = defaults$`output-dir`,
    libraries = libraries,
    merge_distance = as.integer(defaults$`merge-distance`),
    bed_starts_zero_based = parse_bool(defaults$`bed-starts-zero-based`)
  )
}

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required input is missing: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}

run_bedtools <- function(arguments, output_path) {
  bedtools <- Sys.which("bedtools")
  if (!nzchar(bedtools)) {
    stop("bedtools is not available on PATH")
  }
  status <- system2(
    command = bedtools,
    args = arguments,
    stdout = output_path,
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop("bedtools failed with status ", status, ": ", paste(arguments, collapse = " "))
  }
}

build_consensus_peaks <- function(cellranger_root, output_dir, libraries, merge_distance) {
  peak_paths <- vapply(
    libraries,
    function(library_id) {
      require_file(file.path(cellranger_root, library_id, "outs", "peaks.bed"))
    },
    character(1)
  )
  peak_tables <- lapply(
    peak_paths,
    function(path) {
      peaks <- fread(path, header = FALSE, select = 1:3)
      setnames(peaks, c("chrom", "start", "end"))
      peaks
    }
  )
  all_peaks <- rbindlist(peak_tables, use.names = TRUE)
  raw_path <- file.path(output_dir, "all_peaks.raw.bed")
  sorted_path <- file.path(output_dir, "all_peaks.sorted.bed")
  consensus_path <- file.path(output_dir, "consensus_peak_universe.bed")
  fwrite(all_peaks, raw_path, sep = "\t", col.names = FALSE)
  run_bedtools(c("sort", "-i", shQuote(raw_path)), sorted_path)
  run_bedtools(
    c("merge", "-d", as.character(merge_distance), "-i", shQuote(sorted_path)),
    consensus_path
  )
  message("Wrote consensus peak universe to ", normalizePath(consensus_path))
  consensus_path
}

quantify_library <- function(
    library_id,
    cellranger_root,
    consensus_bed_path,
    output_dir,
    bed_starts_zero_based) {
  outs_dir <- file.path(cellranger_root, library_id, "outs")
  fragment_path <- require_file(file.path(outs_dir, "fragments.tsv.gz"))
  require_file(paste0(fragment_path, ".tbi"))
  singlecell_path <- require_file(file.path(outs_dir, "singlecell.csv"))

  message("Reading Cell Ranger barcodes for ", library_id)
  cell_data <- read.csv(singlecell_path, header = TRUE, row.names = 1, check.names = FALSE)
  if (!"is__cell_barcode" %in% colnames(cell_data)) {
    stop(singlecell_path, " does not contain is__cell_barcode")
  }
  filtered_cells <- rownames(cell_data)[cell_data$is__cell_barcode == 1]
  if (length(filtered_cells) == 0L) {
    stop("No Cell Ranger cell barcodes were found for ", library_id)
  }

  consensus_peaks <- fread(
    consensus_bed_path,
    header = FALSE,
    col.names = c("chrom", "start", "end")
  )
  consensus_granges <- makeGRangesFromDataFrame(
    consensus_peaks,
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = bed_starts_zero_based
  )

  message(
    "Quantifying ", length(filtered_cells), " cells across ",
    length(consensus_granges), " consensus peaks for ", library_id
  )
  fragment_object <- CreateFragmentObject(
    path = fragment_path,
    cells = filtered_cells
  )
  counts <- FeatureMatrix(
    fragments = fragment_object,
    features = consensus_granges,
    cells = filtered_cells
  )

  prefix <- file.path(output_dir, paste0(library_id, "_q"))
  writeMM(counts, file = paste0(prefix, "_matrix.mtx"))
  write.table(
    colnames(counts),
    file = paste0(prefix, "_barcodes.tsv"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  fwrite(
    consensus_peaks[, .(chrom, start, end)],
    paste0(prefix, "_peaks.bed"),
    sep = "\t",
    col.names = FALSE
  )
  message("Completed ", library_id, ": ", nrow(counts), " peaks x ", ncol(counts), " cells")
}

main <- function() {
  options(stringsAsFactors = FALSE)
  opt <- parse_args()
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(opt$output_dir, mustWork = TRUE)
  cellranger_root <- normalizePath(opt$cellranger_root, mustWork = TRUE)

  consensus_path <- build_consensus_peaks(
    cellranger_root,
    output_dir,
    opt$libraries,
    opt$merge_distance
  )
  for (library_id in opt$libraries) {
    quantify_library(
      library_id,
      cellranger_root,
      consensus_path,
      output_dir,
      opt$bed_starts_zero_based
    )
  }
  message("Completed consensus-peak requantification for ", length(opt$libraries), " libraries")
}

main()
