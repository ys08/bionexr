# Perform mutation scoring.
# 
# Perform mutation scoring.
# 
# @param x A mutation dataframe.
# @return A dataframe with 3 columns: sample ID, gene name, and score.
# @export
score_mutation <- function(x) {
  ma_file <- paste0(path.package("bionexr"), "/madb.rds")
  if (!file.exists(ma_file)) {
    stop("Please run prepare_ma() first!", call. = FALSE)
  }
  mscore_engine = mutation_assessor(ma_file)
  ret <- mscore_engine(x)
  ret
}

# A simple wrapper for mutation assessor database.
# 
# A simple wrapper for mutation assessor database.
# 
# @param madb_file /inst/madb.rds.
# @return A function for querying mutation assessor database.
# @examples 
# 
# \dontrun{
# madb <- mutation_assessor("inst/madb.rds")
# }
# @export
mutation_assessor <- function(madb_file) {
  stopifnot(file.exists(madb_file))
  
  required_columns <- c("Hugo_Symbol",
                        "Chromosome",
                        "Start_Position",
                        "End_Position",
                        "Reference_Allele",
                        "Tumor_Seq_Allele2",
                        "Tumor_Sample_Barcode")
  function(x) {
    if (ncol(x) != 7) {
      stop("Incorrect mutation data format: ", paste(required_columns, collapse = " "), 
           call. = FALSE)
    }
    
    message("Loading madb database...")
    madb <- readRDS(madb_file)
    search_key <- paste("hg19", x[[2]], x[[3]], x[[5]], x[[6]], sep = "-")
    score <- fast_search_madb(search_key, madb)
    
    ret <- data.frame(sample_id = x[[7]],
                      gene_name = x[[1]],
                      ma_score = score,
                      stringsAsFactors = FALSE)
    ret
  }
}

# Reshape long mutation dataframe to wide mutation dataframe.
# 
# Reshape long mutation dataframe to wide mutation dataframe.
# 
# @param res Long dataframe.
# @param fun Aggregate function.
# @return A wide mutation dataframe.
reshape_long_to_wide <- function(res, fun = max) {
  ret <- reshape2::dcast(res, 
               gene_name ~ sample_id, 
               fun.aggregate = fun, 
               fill = NA_real_)
  row.names(ret) <- ret$gene_name
  ret$gene_name <- NULL
  invisible(ret)
}

# Identify genes with mutations biased towards high functional impact.
# 
# Choose mutated genes above mcutoff and mutation proportion greater than fraction, return the selected genes.
# 
# @param df_mut_score A wide mutation dataframe.
# @param mcutoff A float number between 0 and 6.
# @param fraction A float number between 0 and 1.
# @return A character vector of gene symbols.
choose_mut_genes <- function(df_mut_score, mcutoff, fraction = 0.25) {
  n <- apply(df_mut_score, 1, function(x) sum(!is.na(x)))
  n_p <- apply(df_mut_score, 1, function(x) {
    s <- sum(x >= mcutoff, na.rm = TRUE)
    s
  })
  
  p <- n_p/n
  genes <- names(p)
  genes[p >= fraction]
}