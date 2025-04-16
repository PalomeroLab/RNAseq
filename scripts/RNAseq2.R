# RNA Sequencing Analysis
library('DESeq2')
source("./rnaseq/DESeqDataSetFromFeatureCounts.R")
source("./rnaseq/plots/EnhancedVolcano.R")
source("./rnaseq/plots/pheatmap.R")

dir <- "./data/RNAseq/20240409_Tet2Rhoa-S1P1_RA"
counts_file <- file.path(dir, "counts_noPon.tsv")
# [1] "./data/RNAseq/20240409_Tet2Rhoa-S1P1_RA/counts_noPon.tsv"
design_file <- file.path(dir, "design_noPon.tsv")
# [1] "./data/RNAseq/20240409_Tet2Rhoa-S1P1_RA/design_noPon.tsv"

output_file <- file.path(dir, "dds_noPon.rds")

# Skip this section if we have already generated the DESeq object
dds <- DESeq(DESeqDataSetFromFeatureCounts(counts_file, design_file))
saveRDS(dds, file = output_file)

# Load the DESeq object
dds <- readRDS(output_file)
summary(dds)
colnames(dds)
# [1] "DMSO_1"       "DMSO_2"       "DMSO_3"       "Fingolimod_1" "Fingolimod_2"
# [6] "Fingolimod_3" "Ozanimod_1"   "Ozanimod_2"   "Ozanimod_3"  

resultsNames(dds)
# [1] "Intercept"                    "condition_Fingolimod_vs_DMSO"
# [3] "condition_Ozanimod_vs_DMSO"  

resLFC.fin <- lfcShrink(dds, coef = "condition_Fingolimod_vs_DMSO")
resLFC.oza <- lfcShrink(dds, coef = "condition_Ozanimod_vs_DMSO")

summary(resLFC.fin)
# 
# out of 24765 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1246, 5%
# LFC < 0 (down)     : 936, 3.8%
# outliers [1]       : 0, 0%
# low counts [2]     : 14368, 58%
# (mean count < 28)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# 
# NULL


summary(resLFC.oza)
# 
# out of 24765 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1780, 7.2%
# LFC < 0 (down)     : 1399, 5.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 12051, 49%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# 
# NULL

# resLFC.pon <- lfcShrink(dds, coef = "condition_Ponesimod_vs_DMSO")
# summary(resLFC.pon)


# make a data frame of the geneides and the ld2fc and a third of the pvalue
mcols(resLFC.fin, use.names = TRUE)
mcols(resLFC.oza, use.names = TRUE)

export_DE_stats <- function(obj, prefix, directory) {
  # Check if the directory exists, if not, create it
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  csv_file_path <- file.path(directory, paste0(prefix, "_DE_stats.csv"))
  rnk_file_path <- file.path(directory, paste0(prefix, ".rnk"))
  rnk_noNA_path <- file.path(directory, paste0(prefix, "_noNA.rnk"))

  # Write the full DE stats table as CSV (with rownames as gene IDs)
  write.csv(obj, csv_file_path, row.names = TRUE)
  
  # Create a second data frame containing only gene IDs (rownames) and log2FoldChange
  rnk_data <- data.frame(
    Geneid = rownames(obj),                       # Gene identifiers from rownames
    log2FoldChange = obj$log2FoldChange        # log2FoldChange from the object
  )

  rnk_noNA <- na.omit(rnk_data)  # Remove rows with NA values
  
  # another rnk but with ony significant genes (pval < 0.05)
  rnk_data_sig <- na.omit(rnk_data[obj$padj < 0.05, ])
  
  # Export the .rnk file (tab-separated)
  write.table(rnk_data, rnk_file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(rnk_noNA, rnk_noNA_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(rnk_data_sig, file.path(directory, paste0(prefix, "_sig.rnk")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Return the file paths for reference
  return(list(csv_file = csv_file_path, rnk_file = rnk_file_path))
}

results_dir <- "./results/RNAseq/RUTH/GSEA/individual_analysis"

export_DE_stats(resLFC.fin, "Fingolimod", results_dir)
export_DE_stats(resLFC.oza, "Ozanimod", results_dir)

results_dir <- "./results/RNAseq/RUTH/plots"
# [1] "./results/RNAseq/RUTH/plots"
save_volcano_plot(resLFC.fin, paste0(results_dir, "/volcano_nolabel_Fingolimod.png"), "Fingolimod")
save_volcano_plot(resLFC.oza, paste0(results_dir, "/volcano_nolabel_Ozanimod.png"), "Ozanimod")

plot_and_save_heatmap(dds, "vsd", 5000, results_dir)
plot_and_save_heatmap(dds, assay_type = "ntd", num_genes = 5000, results_dir = results_dir)
plot_and_save_heatmap(dds, assay_type = "rld", num_genes = 5000, results_dir = results_dir)


# If we lower the false discovery rate threshold, we should also inform the
# results() function about it, so that the function can use this threshold for
# the optimal independent filtering that it performs:
res.05 <- results(dds, alpha = 0.05)
mcols(res.05, use.names = TRUE)

# save as csv file with a .rnk extension
write.csv(as.data.frame(res.05), file = "./res_05.csv", row.names = TRUE)



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

# Create subset for DMSO + Fingolimod
dds_fingo <- dds[, dds$condition %in% c("DMSO", "Fingolimod")]
dds_fingo$condition <- factor(dds_fingo$condition, levels = c("DMSO", "Fingolimod"))

# Create VST object for DMSO + Fingolimod comparison
vsd_fingo <- vst(dds_fingo, blind = FALSE)

# Get top variable genes
topVarGenes_fingo <- head(order(rowVars(assay(vsd_fingo)), decreasing = TRUE), 50)
mat_fingo <- assay(vsd_fingo)[topVarGenes_fingo, ]
mat_fingo <- mat_fingo - rowMeans(mat_fingo)

# Create annotation for DMSO + Fingolimod
anno_fingo <- data.frame(condition = colData(vsd_fingo)$condition)
rownames(anno_fingo) <- colnames(vsd_fingo)

# Generate and save heatmap for DMSO + Fingolimod
pheatmap(mat_fingo, 
         annotation_col = anno_fingo, 
         filename = "./heatmap_DMSO_vs_Fingolimod.png", 
         width = 7, 
         height = 7)

# Create subset for DMSO + Ozanimod
dds_ozan <- dds[, dds$condition %in% c("DMSO", "Ozanimod")]
dds_ozan$condition <- factor(dds_ozan$condition, levels = c("DMSO", "Ozanimod"))

# Create VST object for DMSO + Ozanimod comparison
vsd_ozan <- vst(dds_ozan, blind = FALSE)

# Get top variable genes
topVarGenes_ozan <- head(order(rowVars(assay(vsd_ozan)), decreasing = TRUE), 50)
mat_ozan <- assay(vsd_ozan)[topVarGenes_ozan, ]
mat_ozan <- mat_ozan - rowMeans(mat_ozan)

# Create annotation for DMSO + Ozanimod
anno_ozan <- data.frame(condition = colData(vsd_ozan)$condition)
rownames(anno_ozan) <- colnames(vsd_ozan)

# Generate and save heatmap for DMSO + Ozanimod
pheatmap(mat_ozan, 
         annotation_col = anno_ozan, 
         filename = "./heatmap_DMSO_vs_Ozanimod.png", 
         width = 7, 
         height = 7)
