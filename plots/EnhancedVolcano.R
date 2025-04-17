library(EnhancedVolcano)
# Plot the most basic volcano plot
# For the most basic volcano plot, only a single data-frame, data-matrix, or
# tibble of test results is required, containing point labels, log2FC, and
# adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the
# default cut-off for P value is 10e-6.
save_volcano_plot <- function(res_object, path, title) {
  library(EnhancedVolcano)

  p <- EnhancedVolcano(res_object,
    # lab = rownames(res_object),
lab = NA,
    x = "log2FoldChange",
    y = "pvalue",
    title = title,
    pCutoff = 10e-16,
    FCcutoff = 1.5,
    pointSize = 2.0,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1
  )

  ggsave(path, plot = p)
}
