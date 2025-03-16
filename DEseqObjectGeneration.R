#!/usr/bin/env Rscript
source("./DESeqDataSetFromFeatureCounts.R")

# DESeq2 RNA-seq Analysis Script
# This script processes RNA-seq count data using DESeq2

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

# Define command line options
option_list <- list(
  make_option(c("-p", "--path"),
    type = "character", default = NULL,
    help = "Directory containing the counts and design files", metavar = "DIRECTORY"
  ),

  # Parse command line arguments
    opt_parser <- OptionParser(
    option_list = option_list,
    description = "Process RNA-seq count data with DESeq2"
  )
opt <- parse_args(opt_parser)


# Main function
run_deseq_analysis <- function(opt) {
  # Determine file paths based on input parameters
  if (!is.null(opt$directory)) {
    dir_path <- opt$directory

    # Ensure directory path ends with a slash
    if (!endsWith(dir_path, "/")) {
      dir_path <- paste0(dir_path, "/")
    }

    # Set default file paths based on directory
    default_counts <- paste0(dir_path, "counts.tsv")
    default_design <- paste0(dir_path, "design.tsv")
    default_output <- paste0(dir_path, "dds.rds")

    # Check if the specific file format exists (e.g., *_counts.tsv)
    dir_files <- list.files(dir_path)
    counts_pattern <- "_counts\\.tsv$"
    design_patterns <- c("_design\\.tsv$", "_data\\.tsv$", "_metadata\\.tsv$")

    # Find matching count files
    counts_files <- dir_files[grepl(counts_pattern, dir_files)]
    if (length(counts_files) == 1 && is.null(opt$counts)) {
      default_counts <- paste0(dir_path, counts_files)
      message("Found count file: ", default_counts)
    } else if (length(counts_files) > 1 && is.null(opt$counts)) {
      message("Multiple count files found. Using default or specify with --counts:")
      for (i in seq_along(counts_files)) {
        message("  ", i, ": ", counts_files[i])
      }
    }

    # Find matching design files
    for (pattern in design_patterns) {
      design_files <- dir_files[grepl(pattern, dir_files)]
      if (length(design_files) == 1 && is.null(opt$metadata)) {
        default_design <- paste0(dir_path, design_files)
        message("Found design file: ", default_design)
        break
      } else if (length(design_files) > 1 && is.null(opt$metadata)) {
        message("Multiple design files found. Using default or specify with --metadata:")
        for (i in seq_along(design_files)) {
          message("  ", i, ": ", design_files[i])
        }
        break
      }
    }
  } else {
    # If no directory is specified, use current directory
    default_counts <- "counts.tsv"
    default_design <- "design.tsv"
    default_output <- "dds.rds"
  }

  # Determine final file paths with precedence to explicit parameters
  counts_file <- ifelse(!is.null(opt$counts), opt$counts, default_counts)
  design_file <- ifelse(!is.null(opt$metadata), opt$metadata, default_design)
  output_file <- ifelse(!is.null(opt$output), opt$output, default_output)

  # Check if files exist
  if (!file.exists(counts_file)) {
    stop("Counts file not found: ", counts_file)
  }
  if (!file.exists(design_file)) {
    stop("Design file not found: ", design_file)
  }

  # Parse the design formula
  design_formula <- as.formula(opt$formula)

  # Print parameters
  message("\nParameters:")
  message("  Counts file: ", counts_file)
  message("  Design file: ", design_file)
  message("  Output file: ", output_file)
  message("  Design formula: ", as.character(design_formula)[2])
  if (!is.null(opt$reference)) {
    message("  Reference level: ", opt$reference)
  }
  message("")

  # If dry run, stop here
  if (opt$dryrun) {
    message("Dry run completed. Exiting without analysis.")
    return(NULL)
  }

  # Create DESeqDataSet
  dds <- DESeqDataSetFromFeatureCounts(counts_file, design_file, design_formula)

  # Set reference level if specified
  if (!is.null(opt$reference) && "condition" %in% colnames(colData(dds))) {
    if (opt$reference %in% levels(dds$condition)) {
      message("Setting reference level to: ", opt$reference)
      dds$condition <- relevel(dds$condition, ref = opt$reference)
    } else {
      warning(
        "Reference level '", opt$reference, "' not found in condition levels: ",
        paste(levels(dds$condition), collapse = ", ")
      )
    }
  }

  # Run DESeq2 analysis
  message("Running DESeq2 analysis...")
  dds <- DESeq(dds)

  # Save results
  message("Saving DESeqDataSet to: ", output_file)
  saveRDS(dds, file = output_file)

  # Generate and save normalized counts
  norm_counts_file <- sub("\\.rds$", "_normalized_counts.tsv", output_file)
  message("Saving normalized counts to: ", norm_counts_file)
  norm_counts <- counts(dds, normalized = TRUE)
  write.table(norm_counts, file = norm_counts_file, sep = "\t", quote = FALSE)

  # Return the dds object
  return(dds)
}

# Run the analysis if not being sourced
if (!interactive()) {
  # Validate required parameters
  if (is.null(opt$directory) && (is.null(opt$counts) || is.null(opt$metadata))) {
    cat("Error: You must specify either a directory containing counts and design files,\n")
    cat("       or explicitly provide both counts and metadata file paths.\n\n")
    print_help(opt_parser)
    quit(status = 1)
  }

  # Run the analysis
  dds <- run_deseq_analysis(opt)

  if (!is.null(dds)) {
    message("Analysis completed successfully!")
  }
}
