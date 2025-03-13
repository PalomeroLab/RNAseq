# Install and load required packages
install.packages(c("VennDiagram", "gridExtra"))
library(VennDiagram)
library(grid)
library(gridExtra)

# First, let's organize genes by their groupings
# For upregulated genes
# Create sets for the intersections
fin_oza_shared <- intersect(upregulated_Fin, upregulated_Oza)
fin_pon_shared <- intersect(upregulated_Fin, upregulated_Pon)
oza_pon_shared <- intersect(upregulated_Oza, upregulated_Pon)
all_shared_up <- intersect(intersect(upregulated_Fin, upregulated_Oza), upregulated_Pon)
fin_only <- setdiff(upregulated_Fin, union(upregulated_Oza, upregulated_Pon))
oza_only <- setdiff(upregulated_Oza, union(upregulated_Fin, upregulated_Pon))
pon_only <- setdiff(upregulated_Pon, union(upregulated_Fin, upregulated_Oza))
fin_oza_not_pon <- setdiff(fin_oza_shared, upregulated_Pon)
fin_pon_not_oza <- setdiff(fin_pon_shared, upregulated_Oza)
oza_pon_not_fin <- setdiff(oza_pon_shared, upregulated_Fin)

# Function to get gene symbols (assuming rownames are gene IDs)
# This is just an example - modify according to your data structure
get_top_genes <- function(gene_set, n = 5) {
  if(length(gene_set) == 0) return("None")
  if(length(gene_set) <= n) return(paste(gene_set, collapse=", "))
  return(paste(gene_set[1:n], collapse=", "))
}

# Create Venn diagram for upregulated genes
venn_up <- venn.diagram(
  x = list(Fingolimod = upregulated_Fin, 
           Ozanimod = upregulated_Oza, 
           Ponesimod = upregulated_Pon),
  category.names = c("Fingolimod", "Ozanimod", "Ponesimod"),
  filename = NULL,  # Don't save to file yet
  output = FALSE,   # Return the plot object
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

# Create a text grob with the legend
legend_text <- paste(
  "All three drugs: ", get_top_genes(all_shared_up), "\n",
  "Fingolimod & Ozanimod: ", get_top_genes(fin_oza_not_pon), "\n",
  "Fingolimod & Ponesimod: ", get_top_genes(fin_pon_not_oza), "\n",
  "Ozanimod & Ponesimod: ", get_top_genes(oza_pon_not_fin), "\n",
  "Fingolimod only: ", get_top_genes(fin_only), "\n",
  "Ozanimod only: ", get_top_genes(oza_only), "\n",
  "Ponesimod only: ", get_top_genes(pon_only), "\n",
  sep = ""
)

legend_grob <- textGrob(legend_text, 
                       just = "left", 
                       gp = gpar(fontsize = 12))

# Create a layout with the Venn diagram and legend
combined_plot_up <- grid.arrange(
  gTree(children = venn_up), 
  legend_grob, 
  ncol = 1,
  heights = c(3, 1)
)

# Save the combined plot
png("upregulated_venn_with_legend.png", width = 3000, height = 4000, res = 300)
grid.draw(combined_plot_up)
dev.off()

# Repeat the process for downregulated genes
# Create sets for the intersections
fin_oza_shared_down <- intersect(downregulated_Fin, downregulated_Oza)
fin_pon_shared_down <- intersect(downregulated_Fin, downregulated_Pon)
oza_pon_shared_down <- intersect(downregulated_Oza, downregulated_Pon)
all_shared_down <- intersect(intersect(downregulated_Fin, downregulated_Oza), downregulated_Pon)
fin_only_down <- setdiff(downregulated_Fin, union(downregulated_Oza, downregulated_Pon))
oza_only_down <- setdiff(downregulated_Oza, union(downregulated_Fin, downregulated_Pon))
pon_only_down <- setdiff(downregulated_Pon, union(downregulated_Fin, downregulated_Oza))
fin_oza_not_pon_down <- setdiff(fin_oza_shared_down, downregulated_Pon)
fin_pon_not_oza_down <- setdiff(fin_pon_shared_down, downregulated_Oza)
oza_pon_not_fin_down <- setdiff(oza_pon_shared_down, downregulated_Fin)

# Create Venn diagram for downregulated genes
venn_down <- venn.diagram(
  x = list(Fingolimod = downregulated_Fin, 
           Ozanimod = downregulated_Oza, 
           Ponesimod = downregulated_Pon),
  category.names = c("Fingolimod", "Ozanimod", "Ponesimod"),
  filename = NULL,
  output = FALSE,
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

# Create a text grob with the legend for downregulated genes
legend_text_down <- paste(
  "All three drugs: ", get_top_genes(all_shared_down), "\n",
  "Fingolimod & Ozanimod: ", get_top_genes(fin_oza_not_pon_down), "\n",
  "Fingolimod & Ponesimod: ", get_top_genes(fin_pon_not_oza_down), "\n",
  "Ozanimod & Ponesimod: ", get_top_genes(oza_pon_not_fin_down), "\n",
  "Fingolimod only: ", get_top_genes(fin_only_down), "\n",
  "Ozanimod only: ", get_top_genes(oza_only_down), "\n",
  "Ponesimod only: ", get_top_genes(pon_only_down), "\n",
  sep = ""
)

legend_grob_down <- textGrob(legend_text_down, 
                            just = "left", 
                            gp = gpar(fontsize = 12))

# Create a layout with the Venn diagram and legend
combined_plot_down <- grid.arrange(
  gTree(children = venn_down), 
  legend_grob_down, 
  ncol = 1,
  heights = c(3, 1)
)

# Save the combined plot
png("downregulated_venn_with_legend.png", width = 3000, height = 4000, res = 300)
grid.draw(combined_plot_down)
dev.off()
