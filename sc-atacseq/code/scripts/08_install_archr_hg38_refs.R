#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
  cat(
    paste(
      "Usage: 08_install_archr_hg38_refs.R",
      "",
      "Install the hg38 Bioconductor packages required by the ArchR annotation stage into the active R environment.",
      sep = "\n"
    )
  )
  quit(status = 0)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

required_packages <- c(
  "GenomicFeatures",
  "AnnotationDbi",
  "BSgenome.Hsapiens.UCSC.hg38",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) == 0) {
  message("All ArchR hg38 reference packages are already installed.")
  quit(status = 0)
}

message("Installing missing ArchR hg38 reference packages:")
message(paste(" -", missing_packages, collapse = "\n"))
BiocManager::install(missing_packages, ask = FALSE, update = FALSE)
