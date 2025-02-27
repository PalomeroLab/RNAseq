#!/usr/bin/env Rscript
library(DESeq2)
library(apeglm)

# /Users/rdn/GitHub/palomerolab/data/RNAseq/featureCounts
# set the wd to /Users/rdn/GitHub/palomerolab
setwd("/home/ubuntu/20241029_Jurkat_PHF6-PHIP-KO_LQ")
counts_file <- "counts/counts_noclone9.txt"
design_file <- "counts/design.tsv"

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

stopifnot(!is.null(design_data))

stopifnot(all(colnames(count_matrix) == design_data$sample))

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
dds <- DESeq(dds)
res <- results(dds)

# check the resultsNames for the intercepted comparisons
resultsNames(dds)

resLFC_PHF6 <- lfcShrink(dds, coef = "condition_PHF6_vs_CTRL", type = "apeglm")
resLFC_PHIP <- lfcShrink(dds, coef = "condition_PHIP_vs_CTRL", type = "apeglm")

# Extracting transformed values
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
head(assay(vsd), 3)

# Principal component plot of the samples
plotPCA(vsd, intgroup = c("condition"))

# Convert to ggplot by assigning the result to a new var
# Get sample names from dds and add them to pcaData
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

# Get the percentage of variance explained by PC1 and PC2
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Assuming you want to color by the rownames in the dds colData
pcaData$SampleName <- rownames(dds@colData)

# Create the PCA plot, coloring by sample names
ggplot(pcaData, aes(PC1, PC2, color = SampleName, shape = condition)) +
  geom_point(size = 3) +  # Plot the points
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()


# -----------------------------------------------------------------------------
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
# library("RColorBrewer")
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
