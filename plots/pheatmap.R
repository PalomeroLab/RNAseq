library("pheatmap")

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
