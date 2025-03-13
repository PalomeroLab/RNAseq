# Install required packages if you don't have them
install.packages("VennDiagram")
library(VennDiagram)

# Extract significantly upregulated genes (adjust thresholds as needed)
# Using padj < 0.05 and log2FoldChange > 1 as an example
upregulated_Fin <- rownames(resLFC_Fin)[which(resLFC_Fin$padj < 0.05 & resLFC_Fin$log2FoldChange > 1)]
upregulated_Oza <- rownames(resLFC_Oza)[which(resLFC_Oza$padj < 0.05 & resLFC_Oza$log2FoldChange > 1)]
upregulated_Pon <- rownames(resLFC_Pon)[which(resLFC_Pon$padj < 0.05 & resLFC_Pon$log2FoldChange > 1)]

# Extract significantly downregulated genes
downregulated_Fin <- rownames(resLFC_Fin)[which(resLFC_Fin$padj < 0.05 & resLFC_Fin$log2FoldChange < -1)]
downregulated_Oza <- rownames(resLFC_Oza)[which(resLFC_Oza$padj < 0.05 & resLFC_Oza$log2FoldChange < -1)]
downregulated_Pon <- rownames(resLFC_Pon)[which(resLFC_Pon$padj < 0.05 & resLFC_Pon$log2FoldChange < -1)]

# Create Venn diagram for upregulated genes
venn.diagram(
  x = list(Fingolimod = upregulated_Fin, 
           Ozanimod = upregulated_Oza, 
           Ponesimod = upregulated_Pon),
  category.names = c("Fingolimod", "Ozanimod", "Ponesimod"),
  filename = 'upregulated_venn.png',
  output = TRUE,
  imagetype = "png",
  height = 3000, 
  width = 3000, 
  resolution = 300,
  compression = "lzw",
  main = "Upregulated Genes",
  main.cex = 2,
  col = c("#440154", "#21918c", "#fde725"),
  fill = c(alpha("#440154", 0.5), alpha("#21918c", 0.5), alpha("#fde725", 0.5)),
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154", "#21918c", "#fde725"),
  cat.fontface = "bold",
  margin = 0.15
)

# Create Venn diagram for downregulated genes
venn.diagram(
  x = list(Fingolimod = downregulated_Fin, 
           Ozanimod = downregulated_Oza, 
           Ponesimod = downregulated_Pon),
  category.names = c("Fingolimod", "Ozanimod", "Ponesimod"),
  filename = 'downregulated_venn.png',
  output = TRUE,
  imagetype = "png",
  height = 3000, 
  width = 3000, 
  resolution = 300,
  compression = "lzw",
  main = "Downregulated Genes",
  main.cex = 2,
  col = c("#440154", "#21918c", "#fde725"),
  fill = c(alpha("#440154", 0.5), alpha("#21918c", 0.5), alpha("#fde725", 0.5)),
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154", "#21918c", "#fde725"),
  cat.fontface = "bold",
  margin = 0.15
)
