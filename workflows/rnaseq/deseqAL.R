#!/usr/bin/env Rscript
counts_file <- "./data/featureCounts/agx_counts.tsv"
design_file <- "./data/featureCounts/agx_data.tsv"

cat(readLines(design_file), sep = "\n")
# ..Tet2RhoaG17V_24hAGX51_1_sorted_markdup.bam	treated	Treated_1
# ..Tet2RhoaG17V_24hAGX51_2_sorted_markdup.bam	treated	Treated_2
# ..Tet2RhoaG17V_24hAGX51_3_sorted_markdup.bam	treated	Treated_3
# ..Tet2RhoaG17V_24hDMSO_1_sorted_markdup.bam	untreated	Untreated_1
# ..Tet2RhoaG17V_24hDMSO_2_sorted_markdup.bam	untreated	Untreated_2
# ..Tet2RhoaG17V_24hDMSO_3_sorted_markdup.bam	untreated	Untreated_3

dds_al <- DESeq(DESeqDataSetFromFeatureCounts(counts_file, design_file))
dds <- dds_al

resultsNames(dds)
# [1] "Intercept"                      "condition_untreated_vs_treated"

colnames(dds)
# [1] "Treated_1"   "Treated_2"   "Treated_3"   "Untreated_1" "Untreated_2" "Untreated_3"

pca_path <- ("./results/ANOUCH/plots/PCA.png")
save_pca_plot(dds, pca_path)

res <- results(dds, name = "condition_untreated_vs_treated")
resLFC <- lfcShrink(dds, coef = "condition_untreated_vs_treated")
plot_dir <- ("./results/ANOUCH/plots/")
save_volcano_plot(resLFC, paste0(plot_dir, "Volcano.png"), "AGX51 Treated vs Untreated")

plot_and_save_ma(resLFC, plot_dir, "MA", "AGX51 Treated vs Untreated")
