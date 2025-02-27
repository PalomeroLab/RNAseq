#!/usr/bin/env Rscript
library(DESeq2)
library(apeglm)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(dplyr)

setwd("/Users/rdn/GitHub/palomerolab")
counts_file <- "./data/RNAseq/featureCounts/ra_counts.tsv"
design_file <- "./data/RNAseq/featureCounts/ra_data.tsv"

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
design_data$sample
colnames(count_matrix)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = data.frame(
    condition = factor(design_data$condition),
    row.names = design_data$sample
  ),
  design = ~condition
)

# inspect the new dds object
dds

# Run the analysis on the object
# dds <- DESeq(dds)
dds <- DESeq(dds, betaPrior = FALSE)


# check the resultsNames for the intercepted comparisons
resultsNames(dds)

res_fin <- results(dds, name = "condition_Fingolimod_vs_DMSO")
res_oza <- results(dds, name = "condition_Ozanimod_vs_DMSO")
res_pon <- results(dds, name = "condition_Ponesimod_vs_DMSO")
res_fin

resLFC_Fingolimod <- lfcShrink(dds, coef = "condition_Fingolimod_vs_DMSO", type = "apeglm")
resLFC_Ozanimod <- lfcShrink(dds, coef = "condition_Ozanimod_vs_DMSO", type = "apeglm")
resLFC_Ponesimod <- lfcShrink(dds, coef = "condition_Ponesimod_vs_DMSO", type = "apeglm")

# Extracting transformed values
vsd <- vst(dds, blind = FALSE)
# Principal component plot of the samples
# plotPCA(vsd, intgroup = c("condition"))
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

ggsave("./results/RNAseq/RUTH/PCA_plot.png")

colnames(dds)
# Plot the most basic volcano plot
# For the most basic volcano plot, only a single data-frame, data-matrix, or
# tibble of test results is required, containing point labels, log2FC, and
# adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the
# default cut-off for P value is 10e-6.

# 3 volcano plots
EnhancedVolcano(resLFC_Fingolimod,
  lab = rownames(resLFC_Fingolimod),
  x = "log2FoldChange",
  y = "pvalue",
  title = "Fingolimod vs DMSO"
)
ggsave("./results/RNAseq/RUTH/VolcanoFingolimod.png")

EnhancedVolcano(resLFC_Ozanimod,
  lab = rownames(resLFC_Ozanimod),
  x = "log2FoldChange",
  y = "pvalue",
  title = "Ozanimod vs DMSO"
)
ggsave("./results/RNAseq/RUTH/VolcanoOzanimod.png")

EnhancedVolcano(resLFC_Ponesimod,
  lab = rownames(resLFC_Ponesimod),
  x = "log2FoldChange",
  y = "pvalue",
  title = "Ponesimod vs DMSO"
)

ggsave("./results/RNAseq/RUTH/VolcanoPonesimod.png")

# Get the top 10 up and downregulated genes for each treatment
top_genes_list <- list(
  Fingolimod = list(
    up = head(resLFC_Fingolimod[order(resLFC_Fingolimod$log2FoldChange, decreasing = TRUE),], 10),
    down = head(resLFC_Fingolimod[order(resLFC_Fingolimod$log2FoldChange, decreasing = FALSE),], 10)
  ),
  Ozanimod = list(
    up = head(resLFC_Ozanimod[order(resLFC_Ozanimod$log2FoldChange, decreasing = TRUE),], 10),
    down = head(resLFC_Ozanimod[order(resLFC_Ozanimod$log2FoldChange, decreasing = FALSE),], 10)
  ),
  Ponesimod = list(
    up = head(resLFC_Ponesimod[order(resLFC_Ponesimod$log2FoldChange, decreasing = TRUE),], 10),
    down = head(resLFC_Ponesimod[order(resLFC_Ponesimod$log2FoldChange, decreasing = FALSE),], 10)
  )
)

# Add Treatment labels and combine the results into one data frame
top_genes <- bind_rows(
  lapply(names(top_genes_list), function(treatment) {
    treatment_data <- top_genes_list[[treatment]]
    
    # Add Treatment labels and bind rows for up and downregulated genes
    up_df <- as.data.frame(treatment_data$up) %>% mutate(Treatment = paste(treatment, "(Up)"))
    down_df <- as.data.frame(treatment_data$down) %>% mutate(Treatment = paste(treatment, "(Down)"))
    
    bind_rows(up_df, down_df)
  })
)

# Plot
ggplot(top_genes, aes(x = reorder(rownames(top_genes), log2FoldChange), y = log2FoldChange, fill = Treatment)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # To make the plot horizontal
  labs(
    title = "Top 10 Up and Down Regulated Genes for Each Treatment",
    x = "Gene",
    y = "Log2 Fold Change"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("Fingolimod (Up)" = "blue", "Fingolimod (Down)" = "red", 
                              "Ozanimod (Up)" = "green", "Ozanimod (Down)" = "orange", 
                              "Ponesimod (Up)" = "purple", "Ponesimod (Down)" = "yellow"))

ggsave("./results/RNAseq/RUTH/TopGenes.png")

# check S1P1 for each
S1P1 <- counts_data["S1pr1",]

# Extract the LFC for S1pr1 from each DESeq2 result
s1pr1_lfc_fin <- resLFC_Fingolimod["S1pr1", "log2FoldChange"]
s1pr1_lfc_oza <- resLFC_Ozanimod["S1pr1", "log2FoldChange"]
s1pr1_lfc_pon <- resLFC_Ponesimod["S1pr1", "log2FoldChange"]

# Combine the results into a data frame
s1pr1_data <- data.frame(
  Treatment = c("Fingolimod", "Ozanimod", "Ponesimod"),
  log2FoldChange = c(s1pr1_lfc_fin, s1pr1_lfc_oza, s1pr1_lfc_pon)
)

# Plot the LFC for S1pr1 across the treatments
ggplot(s1pr1_data, aes(x = Treatment, y = log2FoldChange, fill = Treatment)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Log2 Fold Change for S1pr1 across Treatments",
    x = "Treatment",
    y = "Log2 Fold Change"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Fingolimod" = "blue", "Ozanimod" = "green", "Ponesimod" = "purple"))


# Extract the LFC and p-values for S1pr1 from each DESeq2 result
s1pr1_lfc_fin <- resLFC_Fingolimod["S1pr1", "log2FoldChange"]
s1pr1_pval_fin <- resLFC_Fingolimod["S1pr1", "pvalue"]

s1pr1_lfc_oza <- resLFC_Ozanimod["S1pr1", "log2FoldChange"]
s1pr1_pval_oza <- resLFC_Ozanimod["S1pr1", "pvalue"]

s1pr1_lfc_pon <- resLFC_Ponesimod["S1pr1", "log2FoldChange"]
s1pr1_pval_pon <- resLFC_Ponesimod["S1pr1", "pvalue"]

# Combine the results into a data frame
s1pr1_data <- data.frame(
  Treatment = c("Fingolimod", "Ozanimod", "Ponesimod"),
  log2FoldChange = c(s1pr1_lfc_fin, s1pr1_lfc_oza, s1pr1_lfc_pon),
  pvalue = c(s1pr1_pval_fin, s1pr1_pval_oza, s1pr1_pval_pon)
)

# Extract the LFC and p-values for S1pr1 from each DESeq2 result
s1pr1_lfc_fin <- resLFC_Fingolimod["S1pr1", "log2FoldChange"]
s1pr1_pval_fin <- resLFC_Fingolimod["S1pr1", "pvalue"]

s1pr1_lfc_oza <- resLFC_Ozanimod["S1pr1", "log2FoldChange"]
s1pr1_pval_oza <- resLFC_Ozanimod["S1pr1", "pvalue"]

s1pr1_lfc_pon <- resLFC_Ponesimod["S1pr1", "log2FoldChange"]
s1pr1_pval_pon <- resLFC_Ponesimod["S1pr1", "pvalue"]

# Combine the results into a data frame
s1pr1_data <- data.frame(
  Treatment = c("Fingolimod", "Ozanimod", "Ponesimod"),
  log2FoldChange = c(s1pr1_lfc_fin, s1pr1_lfc_oza, s1pr1_lfc_pon),
  pvalue = c(s1pr1_pval_fin, s1pr1_pval_oza, s1pr1_pval_pon)
)

# Plot the LFC for S1pr1 across the treatments
ggplot(s1pr1_data, aes(x = Treatment, y = log2FoldChange, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("p=%.5f", pvalue)), vjust = -0.5) +  # Show p-value with 5 decimal places
  labs(
    title = "Log2 Fold Change for S1pr1 across Treatments",
    x = "Treatment",
    y = "Log2 Fold Change"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Fingolimod" = "blue", "Ozanimod" = "green", "Ponesimod" = "purple")) +
  coord_cartesian(ylim = c(min(s1pr1_data$log2FoldChange) - 1, max(s1pr1_data$log2FoldChange) + 1)) +  # Set y-axis limits to center 0
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
ggsave("./results/RNAseq/RUTH/S1pr1_LFC.png")
