library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)

#' Plot PCA and save to file
save_pca_plot <- function(dds, path) {
  # Extract transformed values (`varianceStabilizingTransformation`)
  # It uses the design formula to calculate the within-group variability (if blind=FALSE)
  # or the across-all-samples variability (if blind=TRUE).
  vsd <- vst(dds, blind = FALSE)
  # HACK: This only works for intercepted condition analysis
  pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
  pcaData$sample <- rownames(dds@colData)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # TODO: make the key for sample names match the shapes in the plot
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, shape = condition, color = sample)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  ggsave(path, plot = p)
}

# Plot the most basic volcano plot
# For the most basic volcano plot, only a single data-frame, data-matrix, or
# tibble of test results is required, containing point labels, log2FC, and
# adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the
# default cut-off for P value is 10e-6.
save_volcano_plot <- function(res_object, path, title) {
  library(EnhancedVolcano)

  # TODO: maket this an optional parameter
  # Genes to label explicitly (even if they are not significant)
  # genes_to_label <- c("PHIP", "PHF6")

  # Create a lab vector: label PHIP and PHF6, and also label significant genes
  # labels <- ifelse(rownames(res_object) %in% genes_to_label,
  #   rownames(res_object),
  #   ifelse(res_object$padj < 0.05 & abs(res_object$log2FoldChange) > 1,
  #     rownames(res_object),
  #     ""
  #   )
  # )

  p <- EnhancedVolcano(res_object,
    title = title,
    x = "log2FoldChange",
    y = "pvalue",
    # selectLab = labels,
    lab = rownames(res_object),
  )

  ggsave(path, plot = p)
}


plot_and_save_heatmap <- function(dds, assay_type = c("ntd", "vsd", "rld"), num_genes, results_dir, ...) {
  assay_type <- match.arg(assay_type)

  select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:num_genes]

  df <- data.frame(condition = as.factor(colData(dds)$condition))

  rownames(df) <- colnames(dds)

  # Initialize the transformed_data variable and title
  if (assay_type == "ntd") {
    transformed_data <- normTransform(dds)
    title <- "Normalized counts (without log transformation)"
  } else if (assay_type == "vsd") {
    transformed_data <- vst(dds, blind = FALSE)
    title <- "Variance stabilizing transformation"
  } else if (assay_type == "rld") {
    transformed_data <- rlog(dds, blind = FALSE)
    title <- "Regularized log transformation"
  }

  # Add the number of genes to the title at the end
  title <- paste(title, ", Top", num_genes, "Genes")

  heatmap_data <- assay(transformed_data)
  heatmap_plot <- pheatmap(heatmap_data[select, ],
    # cluster_rows = FALSE,
    # show_rownames = FALSE,
    # cluster_cols = FALSE,
    annotation_col = df,
    silent = TRUE,
    main = title,
    ...
  )

  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  ggsave(paste0(results_dir, "heatmap_", assay_type, num_genes, ".png"),
    plot = heatmap_plot$gtable, width = 8, height = 8
  )
}


plot_and_save_ma <- function(res_object, results_dir, file_prefix, title) {
  png(paste0(results_dir, "/", file_prefix, ".png"), width = 1800, height = 1800, res = 300)
  plotMA(res_object, ylim = c(-2, 2), main = title)
  dev.off()
}
