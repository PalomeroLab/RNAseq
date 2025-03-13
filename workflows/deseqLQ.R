#!/usr/bin/env Rscript
# TODO: maunally check top sig genes in all 3
counts_file <- "./data/LAURA/PHIP/counts_noclone9.txt"
design_file <- "./data/LAURA/PHIP/design_noclone9.csv"
dds_lq <- DESeq(DESeqDataSetFromFeatureCounts(counts_file, design_file))
dds <- dds_lq

cat(readLines(design_file), sep = "\n")
# ..DMSO1_sorted_markdup.bam	DMSO	DMSO_1
# ..DMSO2_sorted_markdup.bam	DMSO	DMSO_2
# ..DMSO3_sorted_markdup.bam	DMSO	DMSO_3
# ..Fingolimod1_sorted_markdup.bam	Fingolimod	Fingolimod_1
# ..Fingolimod2_sorted_markdup.bam	Fingolimod	Fingolimod_2
# ..Fingolimod3_sorted_markdup.bam	Fingolimod	Fingolimod_3
# ..Ozanimod1_sorted_markdup.bam	Ozanimod	Ozanimod_1
# ..Ozanimod2_sorted_markdup.bam	Ozanimod	Ozanimod_2
# ..Ozanimod3_sorted_markdup.bam	Ozanimod	Ozanimod_3
# ..Ponesimod1_sorted_markdup.bam	Ponesimod	Ponesimod_1
# ..Ponesimod2_sorted_markdup.bam	Ponesimod	Ponesimod_2
# ..Ponesimod3_sorted_markdup.bam	Ponesimod	Ponesimod_3

colnames(dds)
# [1] "CTRL_1" "CTRL_2" "CTRL_3" "PHF6_1" "PHF6_2" "PHF6_3" "PHIP_1" "PHIP_2" "PHIP_3"

resultsNames(dds_lq)
# [1] "Intercept"     "condition_PHF6_vs_CTRL" "condition_PHIP_vs_CTRL"

# Using lfcShrink(dds, coef=2)
res_phf6 <- results(dds, name = "condition_PHF6_vs_CTRL") # apeglm
res_phip <- results(dds, name = "condition_PHIP_vs_CTRL") # apeglm

results_dir <- ("results/LAURA/PHIP/no_clone9/MA/")
plot_and_save_ma(resLFC_PHF6_ashr, results_dir, "MA_PHF6_ashr", "PHF6 (ashr)")
plot_and_save_ma(resLFC_PHF6_apeglm, results_dir, "MA_PHF6_apeglm", "PHF6 (apeglm)")
plot_and_save_ma(resLFC_PHF6_normal, results_dir, "MA_PHF6_normal", "PHF6 (normal)")
plot_and_save_ma(resLFC_PHIP_ashr, results_dir, "MA_PHIP_ashr", "PHIP (ashr)")
plot_and_save_ma(resLFC_PHIP_apeglm, results_dir, "MA_PHIP_apeglm", "PHIP (apeglm)")
plot_and_save_ma(resLFC_PHIP_normal, results_dir, "MA_PHIP_normal", "PHIP (normal)")

# Use a function to plot and save the volcano plots automatically
# TODO: make sure PHF6 and PHIP labelled
results_dir <- ("results/LAURA/PHIP/no_clone9/Volcano/labeled/")
plot_and_save_volcano(resLFC_PHF6_ashr,   paste0(results_dir, "/VolcanoPHF6_ashr.png"),   "PHF6 (ashr)")
plot_and_save_volcano(resLFC_PHF6_apeglm, paste0(results_dir, "/VolcanoPHF6_apeglm.png"), "PHF6 (apeglm)")
plot_and_save_volcano(resLFC_PHF6_normal, paste0(results_dir, "/VolcanoPHF6_normal.png"), "PHF6 (normal)")
plot_and_save_volcano(resLFC_PHIP_ashr,   paste0(results_dir, "/VolcanoPHIP_ashr.png"),   "PHIP (ashr)")
plot_and_save_volcano(resLFC_PHIP_apeglm, paste0(results_dir, "/VolcanoPHIP_apeglm.png"), "PHIP (apeglm)")
plot_and_save_volcano(resLFC_PHIP_normal, paste0(results_dir, "/VolcanoPHIP_normal.png"), "PHIP (normal)")


results_dir <- ("results/LAURA/PHIP/no_clone9/heatmap_clustered/")
plot_and_save_heatmap(dds, "ntd", 50, results_dir,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
)

plot_and_save_heatmap(dds, "vsd", 50, results_dir,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
)

plot_and_save_heatmap(dds, "rld", 50, results_dir,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
)
