.libPaths()

counts_file
# [1] "./data/LAURA/PHIP/counts_noclone9.txt"

counts_data <- read.table(
  counts_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = "Geneid"
)

colnames(counts_data[5])


# Extracting transformed values
ntd <- normTransform(dds)
rld <- rlog(dds, blind = FALSE)
vsd <- vst(dds, blind = FALSE)

# Effects of transformations on the variance
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


resLFC_PHIP_ashr <- lfcShrink(dds, coef = "condition_PHIP_vs_CTRL", type = "ashr")
resLFC_PHF6_ashr <- lfcShrink(dds, coef = "condition_PHF6_vs_CTRL", type = "ashr")
resLFC_PHIP_apeglm <- lfcShrink(dds, coef = "condition_PHIP_vs_CTRL", type = "apeglm")
resLFC_PHF6_apeglm <- lfcShrink(dds, coef = "condition_PHF6_vs_CTRL", type = "apeglm")
resLFC_PHIP_normal <- lfcShrink(dds, coef = "condition_PHIP_vs_CTRL", type = "normal")
resLFC_PHF6_normal <- lfcShrink(dds, coef = "condition_PHF6_vs_CTRL", type = "normal")
