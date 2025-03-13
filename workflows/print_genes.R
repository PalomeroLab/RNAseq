# Define a function to pretty-print top N up and down regulated genes and save to a text file
pretty_print_genes <- function(result_list, result_names, results_dir, N) {
  # Create the file name dynamically
  output_file <- file.path(results_dir, paste0("genes", N, ".txt"))

  # Open a connection to the output file
  sink(output_file)

  # Loop through each result dataset
  for (i in 1:length(result_list)) {
    res <- result_list[[i]]
    res_name <- result_names[i]

    original_max_print <- getOption("max.print")
    options(max.print = 10000) # Adjust max.print value as needed


    # Print the name of the variable (dataset)
    cat("\nDataset:", res_name, "\n")

    # Print top N downregulated genes (ordered by increasing log2FoldChange)
    cat(paste("Top", N, "Downregulated Genes:\n"))
    downregulated_genes <- head(res[order(res$log2FoldChange), ], N)
    print(downregulated_genes)

    # Print top N upregulated genes (ordered by decreasing log2FoldChange)
    cat(paste("Top", N, "Upregulated Genes:\n"))
    upregulated_genes <- head(res[order(-res$log2FoldChange), ], N)
    print(upregulated_genes)
  }

  # Close the connection to the output file
  sink()

  options(max.print = original_max_print)

  # Optionally print a message indicating where the file was saved
  cat("Results saved to:", output_file, "\n")
}

# Example of calling the function
result_list <- list(resLFC_PHF6_ashr, resLFC_PHF6_apeglm, resLFC_PHF6_normal, resLFC_PHIP_ashr, resLFC_PHIP_apeglm, resLFC_PHIP_normal)
result_names <- c(
  "resLFC_PHF6_ashr", "resLFC_PHF6_apeglm", "resLFC_PHF6_normal",
  "resLFC_PHIP_ashr", "resLFC_PHIP_apeglm", "resLFC_PHIP_normal"
)

# Define the directory to save the results
results_dir <- ("results/LAURA/PHIP/no_clone9")
# Call the function and save the output to a dynamically named text file with N = 15 top genes
pretty_print_genes(result_list, result_names, results_dir, 50)
