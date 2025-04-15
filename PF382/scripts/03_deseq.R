#!/usr/bin/env Rscript

#' Create a DESeqDataSet from featureCounts output and experimental design data
#' @param counts_file Path to featureCounts output
#' @param design_file Path to two-column design file (no header, any delimiter)
#' @return DESeqDataSet object
#' @import DESeq2
#' @importFrom utils read.table
#' @usage dds <- DESeq(DESeqDataSetFromFeatureCounts(counts_file, design_file))
DESeqDataSetFromFeatureCounts <- function(counts_file, design_file) {
  stopifnot(require(DESeq2))
  stopifnot(all(file.exists(c(counts_file, design_file))))

  tryCatch(
    {
      counts_data <- read.table(
        counts_file,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = "Geneid"
      )

      count_matrix <- as.matrix(counts_data[, 6:ncol(counts_data)])
      mode(count_matrix) <- "integer"

      design_data <- NULL
      for (sep in c("\t", " ", ",")) {
        design_data <- tryCatch(
          read.table(
            design_file,
            sep = sep,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("sample", "condition", "alias")
          ),
          error = function(e) NULL
        )
        if (!is.null(design_data)) {
          if (ncol(design_data) == 2) {
            design_data$alias <- NA
          }
          break
        }
      }

      stopifnot(!is.null(design_data))
      stopifnot(all(colnames(count_matrix) == design_data$sample))

      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = data.frame(
          condition = factor(design_data$condition),
          row.names = design_data$sample
        ),
        design = ~condition
      )

      if (ncol(design_data) == 3) {
        colnames(dds) <- design_data$alias
      }

      mcols(dds)$basepairs <- counts_data[, 5]

      return(dds)
    },
    error = function(e) {
      stop("Error processing files: ", e$message)
    }
  )
}

# Run CLI mode only if this is not sourced
if (!interactive() && identical(environment(), globalenv())) {
  suppressPackageStartupMessages(library(DESeq2))

  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1)

  input_path <- args[1]
  if (!grepl("/$", input_path)) input_path <- paste0(input_path, "/")
  counts_file <- paste0(input_path, "counts.tsv")
  design_file <- paste0(input_path, "design.tsv")
  output_file <- paste0(input_path, "dds.rds")
  norm_counts_file <- paste0(input_path, "dds_normalized_counts.tsv")

  stopifnot(file.exists(counts_file))
  stopifnot(file.exists(design_file))

  message("Creating DESeqDataSet...")
  dds <- DESeqDataSetFromFeatureCounts(counts_file, design_file)

  message("Running DESeq2...")
  dds <- DESeq(dds)

  message("Saving output to: ", output_file)
  saveRDS(dds, file = output_file)

  message("Saving normalized counts to: ", norm_counts_file)
  write.table(counts(dds, normalized = TRUE), file = norm_counts_file, sep = "\t", quote = FALSE)

  message("Done.")
}
