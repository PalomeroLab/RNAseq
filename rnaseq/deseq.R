#!/usr/bin/env Rscript
library(DESeq2)
library(apeglm)

# /Users/rdn/GitHub/palomerolab/data/RNAseq/featureCounts
# set the wd to /Users/rdn/GitHub/palomerolab
setwd("/Users/rdn/GitHub/palomerolab")
counts_file <- "data/RNAseq/featureCounts/lq_counts.tsv"
design_file <- "data/RNAseq/featureCounts/lq_data.tsv"

counts_data <- read.table(
  counts_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = "Geneid"
)

count_matrix <- as.matrix(counts_data[, 6:ncol(counts_data)])
mode(count_matrix) <- "integer"

# Try different separators for design file
design_data <- NULL
for (sep in c("\t", " ", ",")) {
  design_data <- tryCatch(
    read.table(
      design_file,
      sep = sep,
      header = FALSE,
      stringsAsFactors = FALSE,
      col.names = c("sample", "condition")
    ),
    error = function(e) NULL
  )
  if (!is.null(design_data)) break
}

if (is.null(design_data)) {
  stop("Could not parse design file with supported separators")
}

# Validate sample names and create DESeqDataSet
if (!all(colnames(count_matrix) == design_data$sample)) {
  stop("Sample names mismatch between counts and design files")
}
colnames(design_data)
design_data$condition

colnames(count_matrix)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = data.frame(
    condition = factor(design_data$condition),
    row.names = design_data$sample
  ),
  design = ~condition
)
dds

# Here is where you would add additional column metadata
# ....

# Run the analysis on the object
# dds <- DESeq(dds)
dds <- DESeq(dds, betaPrior = FALSE)
res <- results(dds)

# MA-plot
plotMA(res, ylim = c(-2, 2))

# check the resultsNames for the intercepted comparisons
resultsNames(dds)

resLFC_PHF6 <- lfcShrink(dds, coef = "condition_PHF6_vs_CTRL", type = "apeglm")
resLFC_PHIP <- lfcShrink(dds, coef = "condition_PHIP_vs_CTRL", type = "apeglm")

plotMA(resLFC_PHF6, ylim = c(-2, 2))
plotMA(resLFC_PHIP, ylim = c(-2, 2))

# Extracting transformed values
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
head(assay(vsd), 3)

# Principal component plot of the samples
plotPCA(vsd, intgroup = c("condition"))

library(ggplot2)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rowbames(vsd)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


colnames(dds)
# Plot the most basic volcano plot
# For the most basic volcano plot, only a single data-frame, data-matrix, or
# tibble of test results is required, containing point labels, log2FC, and
# adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the
# default cut-off for P value is 10e-6.
library(EnhancedVolcano)

EnhancedVolcano(resLFC_PHIP,
  lab = rownames(resLFC_PHIP),
  x = "log2FoldChange",
  y = "pvalue"
)


EnhancedVolcano(resLFC_PHF6,
  lab = rownames(resLFC_PHF6),
  x = "log2FoldChange",
  y = "pvalue"
)


# Effects of transformations on the variance
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

head(vsdg

