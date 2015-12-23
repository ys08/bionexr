# Convert firehose compressed file to R dataframe.
# 
# Convert firehose compressed file to R dataframe.
# 
# @param file_name firehose tar.gz file.
# @param data_type "stddata_mutation", "stddata_expression", or "analyses_Aggregate_AnalysisFeatures".
# @return An R dataframe object.
# @examples 
# \dontrun{
# d <- dataframe_from_firehose("gdac.broadinstitute.org_HNSC.Mutation_Packager_Calls.Level_3.2015082100.0.0.tar.gz", "stddata_mutation")
# }
# @export
dataframe_from_firehose <- function(file_name, data_type) {
  # check if the compressed file contains MANIFEST.txt
  manifest_pattern <- "MANIFEST\\.txt"
  contents <- untar(file_name, list = TRUE)
  idx <- grepl(manifest_pattern, basename(contents))
  if (sum(idx) != 1) {
    stop("Incorrect firehose file", call. = FALSE)
  }
  # tmp exdir, as short as possible
  tmp_exdir <- paste0("t", Sys.getpid())
  on.exit(unlink(tmp_exdir, recursive = TRUE))
  untar(file_name, exdir = tmp_exdir)
  # exclude manifest
  contents <- contents[!idx]
  # exclude folder end with /
  idx2 <- grepl("/$", contents)
  contents <- contents[!idx2]
  # all files except MANIFEST.txt
  contents <- paste(tmp_exdir, contents, sep = "/")

  switch(data_type, 
         # parse stddata ===================
         # multiple files per tar.gz
         stddata_mutation = {
           df <- dataframe_from_maf(contents)
           invisible(df)
         },
         # one file per tar.gz
         stddata_expression = {
           df <- dataframe_from_rnaseqv2(contents)
           invisible(df)
         },
         #parse analyses data ===================
         analyses_Aggregate_AnalysisFeatures = {
           idx_u <- grepl("\\.samplefeatures", basename(contents))
           if (sum(idx_u == 1)) {
             df <- read.table(contents[idx_u], 
                              header = TRUE, 
                              sep = "\t", 
                              quote = "", 
                              as.is = TRUE, 
                              check.names = FALSE, 
                              fill = TRUE)
             sample_ids <- df[[1]]
             sample_ids <- paste0(sample_ids, "-01")
             sample_idx <- grepl("^TCGA", sample_ids)
             
             df[[1]] <- sample_ids
             df <- df[sample_idx, ]
             
             invisible(df)
           } else {
             stop("Incorrect analyses file", call. = FALSE)
           }
         },
         # default action =================
         {
           stop("Incorrect data type", call. = FALSE)
         })
}

# Create R dataframe from firehose maf file.
# 
# Create R dataframe from firehose maf file.
# 
# Maf file has too many columns, select important columns, keep SNP only and
# return a dataframe.
# 
# @param x A character vector of maf files.
# @return An R dataframe object.
dataframe_from_maf <- function(x) {
  message("Reading and merging ", length(x), " maf files")
  # Hugo_Symbol, Chromosome, Start_position, End_position, 
  # Variant_Type, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode
  ret <- parse_maf(x, col.selected = c(1, 5, 6, 7, 10, 11, 13, 16))
  # only keep SNP type
  ret <- ret[ret[[5]] == "SNP", ]
  ret[5] <- NULL
  ret[[7]] <- substr(ret[[7]], 1, 15)
  # remove duplicated mutation
  dup_idx <- duplicated(ret)
  ret <- ret[!dup_idx, , drop = FALSE]
  ret <- ret[order(ret[[2]], ret[[1]]), , drop = FALSE]
  
  rownames(ret) <- NULL
  invisible(ret)
}

# Create R dataframe from firehose rnaseqv2 file.
# 
# Create R dataframe from firehose rnaseqv2 file.
# 
# @param x Firehose rnaseqv2 file name.
# @return An R dataframe object.
dataframe_from_rnaseqv2 <- function(x) {
  file_size <- get_local_file_size(x)
  file_size <- human_readable(file_size)
  message("Reading ", basename(x), " of size ", file_size, ", and it may take time to finish")
  
  ret <- read.table(x, header = TRUE, sep = "\t", quote = "", as.is = TRUE, fill = TRUE, check.names = FALSE)
  ret <- ret[-1, ]
  
  ret_col <- seq(2, ncol(ret), by = 3)
  gene_id <- gsub("(.+)\\|.+", "\\1", ret[[1]])
  undefined_row <- grepl("\\?", gene_id)
  duplicated_row <- duplicated(gene_id)
  ret_row <- !undefined_row & !duplicated_row
  
  ret <- ret[ret_row, ret_col]
  row.names(ret) <- gene_id[ret_row]
  # adjust sample id
  sample_ids <- colnames(ret)
  colnames(ret) <- substr(sample_ids, 1, 15)
  
  # convert character cell to numeric
  ret[] <- lapply(ret, as.numeric)
  invisible(ret)
}

# Parse firehose maf file.
# 
# Parse firehose maf file.
parse_maf <- function(x, col.selected) {
  stopifnot(!missing(x), !missing(col.selected))
  ret <- lapply(x, function(f) {
    df <- read.table(f, header = T, fill = T, quote = "", as.is = T, sep = "\t")
    df <- df[, col.selected, drop = FALSE]
  })
  
  ret <- Reduce(rbind, ret)
  rownames(ret) <- NULL
  invisible(ret)
}
