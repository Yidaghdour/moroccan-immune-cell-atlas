#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

parse_args <- function(args) {
  parsed <- list(
    peaks_bed = "data/processed/atac_requantified/consensus_peak_universe.bed",
    genome = "hg38",
    output_dir = "results/annotations"
  )

  if ("--help" %in% args || "-h" %in% args) {
    cat(
      paste(
        "Usage: 09_annotate_peaks_archr.R [--peaks-bed PATH] [--genome hg38] [--output-dir DIR]",
        "",
        "Annotate the consensus peak BED with closest genes using ArchR hg38 annotations.",
        "",
        "Options:",
        "  --peaks-bed PATH   Consensus peak BED or BED.gz. Default: data/processed/atac_requantified/consensus_peak_universe.bed",
        "  --genome NAME      ArchR genome preset. Default: hg38",
        "  --output-dir DIR   Output directory. Default: results/annotations",
        sep = "\n"
      )
    )
    quit(status = 0)
  }

  if (length(args) == 0) {
    return(parsed)
  }

  if (length(args) %% 2 != 0) {
    stop("Arguments must be supplied as --key value pairs. Use --help for usage.")
  }

  idx <- seq(1, length(args), by = 2)
  for (i in idx) {
    key <- args[[i]]
    value <- args[[i + 1]]
    if (key == "--peaks-bed") {
      parsed$peaks_bed <- value
    } else if (key == "--genome") {
      parsed$genome <- value
    } else if (key == "--output-dir") {
      parsed$output_dir <- value
    } else {
      stop(sprintf("Unknown argument: %s", key))
    }
  }

  parsed
}

resolve_input_path <- function(path) {
  if (file.exists(path)) {
    return(normalizePath(path, mustWork = TRUE))
  }
  gz_path <- paste0(path, ".gz")
  if (file.exists(gz_path)) {
    return(normalizePath(gz_path, mustWork = TRUE))
  }
  stop(sprintf("Could not find input file %s or %s", path, gz_path))
}

open_text_connection <- function(path, open = "rt") {
  resolved <- resolve_input_path(path)
  if (endsWith(resolved, ".gz")) {
    gzfile(resolved, open = open)
  } else {
    file(resolved, open = open)
  }
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    lines <- c(
      "Missing required R packages for ArchR peak annotation:",
      paste0("  - ", missing),
      "",
      "Install them with:",
      "  make install_archr_refs"
    )
    stop(paste(lines, collapse = "\n"), call. = FALSE)
  }
}

read_peaks <- function(path) {
  con <- open_text_connection(path, open = "rt")
  on.exit(close(con), add = TRUE)
  peaks <- utils::read.table(
    con,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (ncol(peaks) < 3) {
    stop("Peak BED must have at least three columns: chrom, start, end.")
  }
  peaks <- peaks[, 1:3]
  colnames(peaks) <- c("chrom", "start", "end")
  peaks$start <- as.integer(peaks$start)
  peaks$end <- as.integer(peaks$end)
  peaks$peak_id <- paste(peaks$chrom, peaks$start, peaks$end, sep = "-")
  peaks$peak_width <- peaks$end - peaks$start
  peaks$center_1based <- floor((peaks$start + peaks$end) / 2) + 1L
  peaks
}

build_gene_tss <- function(genes_gr) {
  tss_gr <- GenomicRanges::resize(genes_gr, width = 1, fix = "start")
  S4Vectors::mcols(tss_gr)$gene_id <- S4Vectors::mcols(genes_gr)$gene_id
  S4Vectors::mcols(tss_gr)$symbol <- S4Vectors::mcols(genes_gr)$symbol
  S4Vectors::mcols(tss_gr)$strand <- as.character(BiocGenerics::strand(genes_gr))
  tss_gr
}

granges_to_bed_df <- function(gr, name_fun) {
  data.frame(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    start = BiocGenerics::start(gr) - 1L,
    end = BiocGenerics::end(gr),
    name = name_fun(gr),
    score = 0L,
    strand = as.character(BiocGenerics::strand(gr)),
    stringsAsFactors = FALSE
  )
}

write_bed <- function(df, path) {
  ensure_dir(dirname(path))
  con <- if (endsWith(path, ".gz")) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

annotate_supported_peaks <- function(peaks_df, tss_gr) {
  peak_gr <- GenomicRanges::GRanges(
    seqnames = peaks_df$chrom,
    ranges = IRanges::IRanges(start = peaks_df$center_1based, width = 1L)
  )

  hits <- GenomicRanges::distanceToNearest(peak_gr, tss_gr, ignore.strand = TRUE)
  query_idx <- S4Vectors::queryHits(hits)
  subject_idx <- S4Vectors::subjectHits(hits)

  nearest_tss <- tss_gr[subject_idx]
  data.frame(
    peak_id = peaks_df$peak_id[query_idx],
    closest_gene_id = as.character(S4Vectors::mcols(nearest_tss)$gene_id),
    closest_gene_symbol = as.character(S4Vectors::mcols(nearest_tss)$symbol),
    closest_gene_strand = as.character(BiocGenerics::strand(nearest_tss)),
    closest_tss = as.integer(BiocGenerics::start(nearest_tss) - 1L),
    distance_to_tss = abs(peaks_df$center_1based[query_idx] - as.integer(BiocGenerics::start(nearest_tss))),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!identical(args$genome, "hg38")) {
    stop(sprintf("This ArchR annotation stage is configured for hg38 only. Received genome=%s", args$genome))
  }

  require_packages(c(
    "ArchR",
    "GenomicRanges",
    "IRanges",
    "S4Vectors",
    "BiocGenerics",
    "GenomicFeatures",
    "AnnotationDbi",
    "BSgenome.Hsapiens.UCSC.hg38",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db"
  ))

  suppressPackageStartupMessages({
    library(ArchR)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(BiocGenerics)
  })

  peaks <- read_peaks(args$peaks_bed)
  ensure_dir(args$output_dir)

  genome_annotation <- ArchR::createGenomeAnnotation(genome = args$genome)
  gene_annotation <- ArchR::createGeneAnnotation(genome = args$genome)

  genes_gr <- gene_annotation$genes
  gene_tss_gr <- build_gene_tss(genes_gr)

  supported_seqlevels <- intersect(
    unique(peaks$chrom),
    intersect(
      unique(as.character(GenomicRanges::seqnames(genome_annotation$chromSizes))),
      unique(as.character(GenomicRanges::seqnames(gene_tss_gr)))
    )
  )

  supported_mask <- peaks$chrom %in% supported_seqlevels
  unsupported_peaks <- peaks[!supported_mask, , drop = FALSE]
  supported_peaks <- peaks[supported_mask, , drop = FALSE]
  supported_tss_gr <- gene_tss_gr[as.character(GenomicRanges::seqnames(gene_tss_gr)) %in% supported_seqlevels]

  annotation_df <- data.frame(
    chrom = peaks$chrom,
    start = peaks$start,
    end = peaks$end,
    peak_id = peaks$peak_id,
    peak_center = peaks$center_1based - 1L,
    peak_width = peaks$peak_width,
    closest_gene_id = NA_character_,
    closest_gene_symbol = NA_character_,
    closest_gene_strand = NA_character_,
    closest_tss = NA_integer_,
    distance_to_tss = NA_integer_,
    stringsAsFactors = FALSE
  )

  if (nrow(supported_peaks) > 0) {
    supported_ann <- annotate_supported_peaks(supported_peaks, supported_tss_gr)
    match_idx <- match(supported_ann$peak_id, annotation_df$peak_id)
    annotation_df$closest_gene_id[match_idx] <- supported_ann$closest_gene_id
    annotation_df$closest_gene_symbol[match_idx] <- supported_ann$closest_gene_symbol
    annotation_df$closest_gene_strand[match_idx] <- supported_ann$closest_gene_strand
    annotation_df$closest_tss[match_idx] <- supported_ann$closest_tss
    annotation_df$distance_to_tss[match_idx] <- supported_ann$distance_to_tss
  }

  gene_bed <- granges_to_bed_df(
    genes_gr,
    name_fun = function(gr) {
      paste(S4Vectors::mcols(gr)$symbol, S4Vectors::mcols(gr)$gene_id, sep = "|")
    }
  )
  tss_bed <- granges_to_bed_df(
    gene_tss_gr,
    name_fun = function(gr) {
      paste(S4Vectors::mcols(gr)$symbol, S4Vectors::mcols(gr)$gene_id, sep = "|")
    }
  )

  annotation_path <- file.path(args$output_dir, "consensus_peak_annotations.tsv.gz")
  gene_bed_path <- file.path(args$output_dir, "reference_genes.bed.gz")
  tss_bed_path <- file.path(args$output_dir, "reference_tss.bed.gz")

  con <- gzfile(annotation_path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(annotation_df, file = con, sep = "\t", quote = FALSE, row.names = FALSE)

  write_bed(gene_bed, gene_bed_path)
  write_bed(tss_bed, tss_bed_path)

  if (nrow(unsupported_peaks) > 0) {
    unsupported_counts <- sort(table(unsupported_peaks$chrom), decreasing = TRUE)
    warning(
      paste0(
        "Peaks on unsupported seqlevels were left unannotated: ",
        paste(sprintf("%s=%s", names(unsupported_counts), unsupported_counts), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  message(sprintf("Wrote peak annotations to %s", normalizePath(annotation_path, mustWork = FALSE)))
  message(sprintf("Wrote gene BED to %s", normalizePath(gene_bed_path, mustWork = FALSE)))
  message(sprintf("Wrote TSS BED to %s", normalizePath(tss_bed_path, mustWork = FALSE)))
  message(sprintf("Annotated %s peaks (%s supported, %s unsupported).", nrow(peaks), nrow(supported_peaks), nrow(unsupported_peaks)))
}

main()
