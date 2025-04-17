#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)

input_path <- args[1]
if (!grepl("/$", input_path)) input_path <- paste0(input_path, "/")
counts_file <- paste0(input_path, "counts.tsv")
design_file <- paste0(input_path, "design.tsv")
output_file <- paste0(input_path, "dds.rds")
norm_counts_file <- paste0(input_path, "dds_normalized_counts.tsv")

script_dir <- dirname(normalizePath(sys.frame(1)$ofile))
source(file.path(script_dir, "DESeqDataSetFromFeatureCounts.R"))

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
