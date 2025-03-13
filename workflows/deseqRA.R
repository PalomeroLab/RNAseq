
colnames(dds)
#  [1] "DMSO_1"       "DMSO_2"       "DMSO_3"       "Fingolimod_1" "Fingolimod_2" "Fingolimod_3"
#  [7] "Ozanimod_1"   "Ozanimod_2"   "Ozanimod_3"   "Ponesimod_1"  "Ponesimod_2"  "Ponesimod_3"

resultsNames(dds)
# "condition_Fingolimod_vs_DMSO" "condition_Ozanimod_vs_DMSO"   "condition_Ponesimod_vs_DMSO"

res_fin <- results(dds, name = "condition_Fingolimod_vs_DMSO")
res_oza <- results(dds, name = "condition_Ozanimod_vs_DMSO")
res_pon <- results(dds, name = "condition_Ponesimod_vs_DMSO")

resLFC_Fin <- lfcShrink(dds, coef = "condition_Fingolimod_vs_DMSO")
resLFC_Oza <- lfcShrink(dds, coef = "condition_Ozanimod_vs_DMSO")
resLFC_Pon <- lfcShrink(dds, coef = "condition_Ponesimod_vs_DMSO")

pca_path <- ("./results/RUTH/plots/PCA.png")
save_pca_plot(dds, pca_path)

plot_dir <- ("./results/RUTH/plots/")
save_volcano_plot(resLFC_Fin, paste0(plot_dir, "VolcanoFingolimod.png"), "Fingolimod vs DMSO")
save_volcano_plot(resLFC_Oza, paste0(plot_dir, "VolcanoOzanimod.png"),   "Ozanimod vs DMSO")
save_volcano_plot(resLFC_Pon, paste0(plot_dir, "VolcanoPonesimod.png"),  "Ponesimod vs DMSO")
# TODO: function to stitch 3 plots


# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
# Check the structure of df to confirm it has the 'condition' column


select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:150]
df <- data.frame(condition = as.factor(colData(dds)$condition))
rownames(df) <- colnames(dds)

# Capture the pheatmap plot for ntd
heatmap_ntd <- pheatmap(
  assay(ntd)[select, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  annotation_col = df,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  silent = TRUE # This prevents pheatmap from showing the plot immediately
)

# Save the pheatmap for ntd using ggsave
ggsave("./results/RNAseq/RUTH/heatmap_normalized50.png", plot = heatmap_ntd$gtable, width = 8, height = 8)


# Capture the pheatmap plot for vsd
heatmap_vsd <- pheatmap(
  assay(vsd)[select, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  annotation_col = df,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  silent = TRUE # This prevents pheatmap from showing the plot immediately
)

# Save the pheatmap for vsd using ggsave
ggsave("./results/RNAseq/RUTH/heatmap_stabalized150.png", plot = heatmap_vsd$gtable, width = 8, height = 8)



# Get the top 10 up and downregulated genes for each treatment
top_genes_list <- list(
  Fingolimod = list(
    up = head(resLFC_Fingolimod[order(resLFC_Fingolimod$log2FoldChange, decreasing = TRUE), ], 10),
    down = head(resLFC_Fingolimod[order(resLFC_Fingolimod$log2FoldChange, decreasing = FALSE), ], 10)
  ),
  Ozanimod = list(
    up = head(resLFC_Ozanimod[order(resLFC_Ozanimod$log2FoldChange, decreasing = TRUE), ], 10),
    down = head(resLFC_Ozanimod[order(resLFC_Ozanimod$log2FoldChange, decreasing = FALSE), ], 10)
  ),
  Ponesimod = list(
    up = head(resLFC_Ponesimod[order(resLFC_Ponesimod$log2FoldChange, decreasing = TRUE), ], 10),
    down = head(resLFC_Ponesimod[order(resLFC_Ponesimod$log2FoldChange, decreasing = FALSE), ], 10)
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
  coord_flip() + # To make the plot horizontal
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
  scale_fill_manual(values = c(
    "Fingolimod (Up)" = "blue", "Fingolimod (Down)" = "red",
    "Ozanimod (Up)" = "green", "Ozanimod (Down)" = "orange",
    "Ponesimod (Up)" = "purple", "Ponesimod (Down)" = "yellow"
  ))

ggsave("./results/RNAseq/RUTH/TopGenes.png")

# check S1P1 for each
S1P1 <- counts_data["S1pr1", ]

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
  geom_text(aes(label = sprintf("p=%.5f", pvalue)), vjust = -0.5) + # Show p-value with 5 decimal places
  labs(
    title = "Log2 Fold Change for S1pr1 across Treatments",
    x = "Treatment",
    y = "Log2 Fold Change"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Fingolimod" = "blue", "Ozanimod" = "green", "Ponesimod" = "purple")) +
  coord_cartesian(ylim = c(min(s1pr1_data$log2FoldChange) - 1, max(s1pr1_data$log2FoldChange) + 1)) + # Set y-axis limits to center 0
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
ggsave("./results/RNAseq/RUTH/S1pr1_LFC.png")
