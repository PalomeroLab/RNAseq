setwd("/Users/rdn/GitHub/palomerolab/workflows/rnaseq/")
source(".Rprofile")

# Combine the function calls for DESeq our wrapper function
# note: betaPrior = FALSE is default
dds <- DESeq(DESeqDataSetFromFeatureCounts(counts_file, design_file))

resultsNames(dds_default)
# [1] "Intercept"     "condition_PHF6_vs_CTRL" "condition_PHIP_vs_CTRL"

resultsNames(dds)
# [1] "Intercept"     "condition_PHF6_vs_CTRL" "condition_PHIP_vs_CTRL"

resultsNames(ddsB)
# [1] "Intercept"     "conditionCTRL" "conditionPHF6" "conditionPHIP"

# Using  lfcShrink(dds, coef=2)
res_phf6 <- results(dds, coef=2, name = "condition_PHF6_vs_CTRL")
res_phip <- results(dds, coef=2, name = "condition_PHIP_vs_CTRL")

# MA-plot
plotMA(res, ylim = c(-2, 2))
