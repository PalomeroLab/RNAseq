#' Create a DESeqDataSet from featureCounts output and experimental design data
#' @param counts_file Path to featureCounts output
#' @param design_file Path to two-column design file (no header, any delimiter)
#' @return DESeqDataSet object
#' @import DESeq2
#' @importFrom utils read.table
#' @usage dds <- DESeq(DESeqDataSetFromFeatureCounts(counts_file, design_file))
DESeqDataSetFromFeatureCounts <- function(counts_file, design_file) {
  stopifnot(require(DESeq2))
  stopifnot(all(file.exists(c(counts_file, design_file))))

  tryCatch(
    {
      counts_data <- read.table(
        counts_file,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = "Geneid"
      )

      # Extract the count matrix (assuming the first 5 columns are metadata)
      count_matrix <- as.matrix(counts_data[, 6:ncol(counts_data)])
      mode(count_matrix) <- "integer"

      # Try different separators for design file
      design_data <- NULL
      for (sep in c("\t", " ", ",")) {
        design_data <- tryCatch(
          read.table(
            design_file,
            sep = sep,
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = c("sample", "condition", "alias")
          ),
          error = function(e) NULL
        )

        if (!is.null(design_data)) {
          # If the file has only two columns,
          # add the alias column with NA
          if (ncol(design_data) == 2) {
            design_data$alias <- NA
          }
          break
        }
      }

      # If no design data could be read, stop the function with an error message
      stopifnot(!is.null(design_data))

      # Validate that the sample names in counts and design files match
      # DESeq2 requires them to be in the same order as they appear in featureCounts output
      stopifnot(all(colnames(count_matrix) == design_data$sample))

      # Create DESeqDataSet from the count matrix and experimental design
      dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = data.frame(
          condition = factor(design_data$condition),
          row.names = design_data$sample
        ),
        design = ~condition
      )

      # if there were aliases, replace the sample names with them
      if (ncol(design_data) == 3) {
        colnames(dds) <- design_data$alias
      }

      # append a column to dds, mcols(object)$basepairs to produce FPKM values
      # colnames(counts_data[5]) is the basepaits
      mcols(dds)$basepairs <- counts_data[, 5]

      return(dds)
    },
    error = function(e) {
      stop("Error processing files: ", e$message)
    }
  )
}
