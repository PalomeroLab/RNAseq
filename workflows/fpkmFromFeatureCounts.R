#!/usr/bin/env Rscript

#' countToFPKM: Convert Counts to Fragments per Kilobase of Transcript per Million (FPKM)
library(ComplexHeatmap)
library(countToFPKM)

# args <- commandArgs(trailingOnly = TRUE)
# counts_file <- args[1]
# stopifnot(file.exists(counts_file))

# inputs# feature counts with m(cols)$basepairs filled and set to feature length
fpkm_data <- DESeq2::fpkm(dds)
png("results/LAURA/PHIP/no_clone9/fpkmheatmap.png", width = 1800, height = 1800, res = 300 )
fpkmheatmap(fpkm_data, topvar=30, showfeaturenames=TRUE, return_log = TRUE)
dev.off()
