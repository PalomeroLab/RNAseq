# RNA Sequencing Analysis
renv::use(lockfile = "./rnaseq/renv.lock")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
library(DESeq2)
library(apeglm)
source("./rnaseq/R/DESeqDataSetFromFeatureCounts.R")

dir <- "./data/RNAseq/20240409_Tet2Rhoa-S1P1_RA"
counts_file <- file.path(dir, "counts_noPon.tsv")
# [1] "./data/RNAseq/20240409_Tet2Rhoa-S1P1_RA/counts_noPon.tsv"
design_file <- file.path(dir, "design_noPon.tsv")
# [1] "./data/RNAseq/20240409_Tet2Rhoa-S1P1_RA/design_noPon.tsv"

counts_file <- file.path(dir, "counts.tsv")
design_file <- file.path(dir, "design.tsv")

output_file <- file.path(dir, "dds_noPon.rds")

# Skip this section if we have already generated the DESeq object
dds <- DESeqDataSetFromFeatureCounts(counts_file, design_file)
dds <- DESeq(DESeqDataSetFromFeatureCounts(counts_file, design_file))
saveRDS(dds, file = output_file)

# Load the DESeq object
dds <- readRDS(output_file)
summary(dds)
colnames(dds)
#  [1] "DMSO_1"       "DMSO_2"       "DMSO_3"       "Fingolimod_1" "Fingolimod_2"
#  [6] "Fingolimod_3" "Ozanimod_1"   "Ozanimod_2"   "Ozanimod_3"   "Ponesimod_1"
# [11] "Ponesimod_2"  "Ponesimod_3,

resultsNames(dds)


resLFC.fin <- lfcShrink(dds, coef = "condition_Fingolimod_vs_DMSO")
resLFC.oza <- lfcShrink(dds, coef = "condition_Ozanimod_vs_DMSO")

summary(resLFC.fin)
# out of 25577 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1306, 5.1%
# LFC < 0 (down)     : 1042, 4.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 14402, 56%
# (mean count < 19)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

summary(resLFC.oza)
# out of 25577 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1833, 7.2%
# LFC < 0 (down)     : 1464, 5.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 12481, 49%
# (mean count < 7)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

resLFC.pon <- lfcShrink(dds, coef = "condition_Ponesimod_vs_DMSO")
summary(resLFC.pon)
# out of 25577 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3880, 15%
# LFC < 0 (down)     : 3931, 15%
# outliers [1]       : 0, 0%
# low counts [2]     : 8641, 34%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# If we lower the false discovery rate threshold, we should also inform the
# results() function about it, so that the function can use this threshold for
# the optimal independent filtering that it performs:
res.05 <- results(dds, alpha = 0.05)
mcols(res.05, use.names = TRUE)
# DataFrame with 6 rows and 2 columns
#                        type            description
#                 <character>            <character>
# baseMean       intermediate mean of normalized c..
# log2FoldChange      results log2 fold change (ML..
# lfcSE               results standard error: cond..
# stat                results Wald statistic: cond..
# pvalue              results Wald test p-value: c..
# padj                results   BH adjusted p-values

table(res.05$padj < 0.05)
#
# FALSE  TRUE
#  8131  6885

resLFC1 <- results(dds, lfcThreshold = 1)
table(resLFC1$padj < 0.05)
#
# FALSE  TRUE
# 16188   268

# Gene clustering
BiocManager::install("pheatmap", package="binary")
library("genefilter")
library("pheatmap")
vsd <- vst(dds, blind = FALSE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
anno <- data.frame(condition = colData(vsd)$condition)
rownames(anno) <- colnames(vsd)
pheatmap(mat, annotation_col = anno, filename = "./heatmap2_100.png", width = 7, height = 7)

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
anno <- data.frame(condition = colData(vsd)$condition)
rownames(anno) <- colnames(vsdl
pheatmap(mat, annotation_col = anno, filename = "./heatmap50.png", width = 10, height = 10)


# inspect S1P1 in each table
res.05[res.05$gene == "S1pr1", ]
# print the values
resLFC.fin[resLFC.fin$gene == "S1pr1", ]
print(resLFC.oza[resLFC.oza$gene == "S1pr1", ])

# Extract rows where gene == "S1pr1" in the res.05 table
s1pr1_res_05 <- res.05[res.05$gene == "S1pr1", ]
print(s1pr1_res_05)

# Extract rows where gene == "S1pr1" in the resLFC.fin table
s1pr1_resLFC_fin <- resLFC.fin[resLFC.fin$gene == "S1pr1", ]
print(s1pr1_resLFC_fin)

# Extract rows where gene == "S1pr1" in the resLFC.oza table
s1pr1_resLFC_oza <- resLFC.oza[resLFC.oza$gene == "S1pr1", ]
print(s1pr1_resLFC_oza)


# check for al genes starting wit 's1p'
for (gene in rownames(res.05)) {
  if (startsWith(gene, "S1p")) {
    print(res.05[res.05$gene == gene, ])
  }
}

# Loop through the rownames of res.05 to check for genes starting with 'S1p'
for (gene in rownames(res.05)) {
  if (startsWith(gene, "S1p")) {
    # Print the gene name
    print(gene)
    # Print the rows corresponding to that gene
    print(res.05[res.05$gene == gene, ])
  }
}

# Check the column names to ensure the 'gene' column exists
print(colnames(res.05))

# Loop through the rownames of res.05 to check for genes starting with 'S1p'
for (gene in rownames(res.05)) {
  if (startsWith(gene, "S1p")) {
    # Check if the 'gene' column contains this gene
    matching_rows <- res.05[res.05$gene == gene, ]
    
    # Print the gene name if it has matching rows
    if (nrow(matching_rows) > 0) {
      print(paste("Gene:", gene))
      print(matching_rows)
    }
  }
}

# Check the structure of the data frame to confirm column names
str(res.05)

# Print the first few rows to check how the data looks
head(res.05)

# Check if the column 'gene' exists
if ("gene" %in% colnames(res.05)) {
  cat("'gene' column exists in res.05\n")
} else {
  cat("'gene' column does not exist in res.05\n")
}

# Check for the row names of the data frame (if you want to match row names with gene)
cat("Row names of res.05: \n")
head(rownames(res.05))

# Loop through the gene names in row names and check for matches
for (gene in rownames(res.05)) {
  if (startsWith(gene, "S1p")) {
    # Find rows where gene == the row name
    matching_rows <- res.05[res.05$gene == gene, ]

    # Only print if matching rows exist
    if (nrow(matching_rows) > 0) {
      print(paste("Gene found:", gene))
      print(matching_rows)
    } else {
      print(paste("No matching rows found for gene:", gene))
    }
  }
}

