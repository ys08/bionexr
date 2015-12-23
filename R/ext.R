# Build a query database for mutation assessor dataset.
# 
# Build a named vector for the mutation assessor dataset. Use data.table
# package to read the large file.
# 
# @param ma_file Mutation assessor dataset.
# @param ... parameters for data.table::fread.
# @return A named vector.
# @examples 
# 
# \dontrun{
# madb <- build_madb(ma_file)
# }
# build_madb <- function(ma_file, ...) {
#   if(!requireNamespace("data.table", quietly = TRUE)) {
#     stop("This function needs the data.table library installed, albeit it cannot be found. Check out the installation instructions!\n", 
#          call. = FALSE)
#   }
#   
#   dt <- data.table::fread(ma_file, ...)
#   ret <- dt$V2
#   names(ret) <- dt$V1
#   
#   stopifnot(is.numeric(ret))
#   
#   invisible(ret)
# }

search_madb <- function(ids, madb = stop("'madb' must be specified")) {
  message("Scoring missense mutation by Mutation Assessor")
  message("It will take about", round(length(ids) * 0.63 / 3600, 1), "hours to finish querying ", length(ids), " mutations")
  
  ret <- vector(mode = "numeric", length = length(ids))
  pb <- txtProgressBar(1, length(ids), style=3)
  for (i in seq_along(ids)) {
    ret[i] <- tryCatch(madb[[ids[i]]], error = function(e) 0)
    setTxtProgressBar(pb, i)
  }
  ret
}

# Query mutation assessor database.
# 
# This is the single thread query function, ~0.37 seconds per query using one 2.7Ghz CPU.
# 
# @param ids query IDs.
# @param madb A named vector of mutation assessor dataset.
# @return The functional impact scores of \code{ids}.
# @export
fast_search_madb_single_thread <- function(ids, madb = stop("'madb' must be specified")) {
  message("Scoring missense mutation by Mutation Assessor")
  message("It will take about ", round(length(ids) * 0.37 / 3600, 1), " hours to finish querying ", length(ids), " mutations")
  
  ret <- vector(mode = "numeric", length = length(ids))
  idx <- ids %in% names(madb)
  pb <- txtProgressBar(1, length(ids), style = 3)
  
  for (i in seq_along(ids)) {
    if (idx[i]) {
      ret[i] <- madb[[ids[i]]]
    } else {
      ret[i] <- 0
    }
    setTxtProgressBar(pb, i)
  }
  ret
}
# Query mutation assessor database.
# 
# This is the multi-thread query function, ~0.37 seconds per query using one 
# 2.7Ghz CPU.
# 
# This function chooses automaticly to run in parallel mode if the host 
# computer's physical memory size is larger than 8GB. It uses "doParallel" 
# package to paralize the query process. To many threads will run into memory 
# issues because "doParallel" copys the madb object to every thread, and madb
# object will eat all of the memmory.
# @param ids query IDs.
# @param madb A named vector of mutation assessor dataset.
# @return The functional impact scores of \code{ids}.
# @export
fast_search_madb <- function(ids, madb = stop("'madb' must be specified")) {
  ret <- vector(mode = "numeric", length = length(ids))
  idx <- ids %in% names(madb)
  N <- length(ids)
  
  if(requireNamespace("doParallel", quietly = TRUE) && memory.limit() >= 80000) {
    n_cores <- 2
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    log_file <- "mutation_scoring_progress_multi.txt"
    write(as.character(Sys.time()), log_file)
    
    message("Scoring missense mutation by Mutation Assessor using ", n_cores, " cores")
    message("It will take about ", round(length(ids) * 0.37 / (3600 * n_cores), 1), " hours to finish querying ", length(ids), " mutations")
    
    ret <- foreach::`%dopar%`(foreach::foreach(i = seq_len(N), .combine = c), {
      sink(log_file, append = TRUE)
      cat(paste0(i, "/", N, "\n"))
      
      res <- 0
      if (idx[i]) {
        res <- madb[[ids[i]]]
      }
      sink()
      res
    })
    
    parallel::stopCluster(cl)
    write(as.character(Sys.time()), log_file, append = TRUE)
    
    return(ret)
  } else {
    message("Scoring missense mutation by Mutation Assessor")
    message("It will take about ", round(length(ids) * 0.37 / 3600, 1), " hours to finish querying ", length(ids), " mutations")

    pb <- txtProgressBar(1, length(ids), style=3)
    
    for (i in seq_along(ids)) {
      if (idx[i]) {
        ret[i] <- madb[[ids[i]]]
      } else {
        ret[i] <- 0
      }
      setTxtProgressBar(pb, i)
    }
    return(ret)
  }
}

# Preprocess user-provided data for PPI-based analysis.
# 
# Preprocess user-provided data for PPI-based analysis.
# 
# Preprocess user-provided data for PPI-based analysis. \code{df_exp}: column 
# names -- sample id; row names -- gene symbols. Write exluded data to 
# "user_mutated_gene_exclude_ppi.txt" and 
# "user_expressed_gene_exclude_ppi.txt".
#
# @param df_mut user-provided mutation data.
# @param df_exp user-provided expression data.
# @param sid_pattern A string regex.
# @return A list, first element is mutation data, second element is expression
#   data, third element is case/control sample index in expression data.
# @examples 
# 
# \dontrun{
# res <- preprocess_user_data_ppi(hnsc_mut, hnsc_exp, c("-01$", "-11$"))
# }
# @export
preprocess_user_data_ppi <- function(df_mut, df_exp, sid_pattern) {
  ret <- vector(mode = "list", length = 3)
  # check user's data column name
  mut_colnames <- colnames(df_mut)
  if (length(mut_colnames) != 7|
      mut_colnames[1] != "Hugo_Symbol"|
      mut_colnames[2] != "Chromosome"|
      mut_colnames[3] != "Start_position"|
      mut_colnames[4] != "End_position"|
      mut_colnames[5] != "Reference_Allele"|
      mut_colnames[6] != "Tumor_Seq_Allele2"|
      mut_colnames[7] != "Tumor_Sample_Barcode") {
    stop("Your gene mutation data is invalid, please check again!", call. = FALSE)
  }
  
  mutated_symbols <- df_mut[[1]]
  expressed_symbols <- rownames(df_exp)
  mut_idx <- mutated_symbols %in% valid_ppi_symbols()
  exp_idx <- expressed_symbols %in% valid_ppi_symbols()
  
  # a few genes are not analyzed and are logged to files.
  write.table(df_mut[!mut_idx, ], file = "user_mutated_gene_exclude_ppi.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(df_exp[!exp_idx, ], file = "user_expressed_gene_exclude_ppi.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  # drop some user's data
  mut_included <- df_mut[mut_idx, ]
  exp_included <- df_exp[exp_idx, ]
  
  case_pattern <- sid_pattern[1]
  control_pattern <- sid_pattern[2]
  
  mut_sid <- unique(mut_included[[7]])
  mut_case <- grepl(case_pattern, mut_sid)
  if (!all(mut_case)) {
    stop("Your mutation data contains invalid sample IDs, please check again!", call. = FALSE)
  }
  
  exp_sid <- colnames(exp_included)
  exp_case <- grepl(case_pattern, exp_sid)
  exp_control <- grepl(control_pattern, exp_sid)
  if (!any(exp_case) | !any(exp_control)) {
    stop("Your expression data does not contain case/control data, please check again!", call. = FALSE)
  }
  
  if (sum(exp_case | exp_control) != length(exp_sid)) {
    stop("Your expression data contains invalid sample IDs, please check again!", call. = FALSE)
  }
  
  if (!setequal(mut_sid, exp_sid[exp_case])) {
    warning("Your mutation data's sample IDs probably does not match your expression data's sample IDs, please take carefull note of the analysis result!", call. = FALSE)
  }
  
  exp_case_sid_sorted <- sort(exp_sid[exp_case])
  exp_control_sid_sorted <- sort(exp_sid[exp_control])
  ret_exp_case_idx <- match(exp_case_sid_sorted, exp_sid)
  ret_exp_control_idx <- match(exp_control_sid_sorted, exp_sid)
  
  # return
  ret[[1]] <- mut_included
  ret[[2]] <- exp_included
  ret[[3]] <- list(case = ret_exp_case_idx, control = ret_exp_control_idx)
  ret
}
# Preprocess user-provided data for pathway-based analysis.
# 
# Preprocess user-provided data for pathway-based analysis.
# 
# Preprocess user-provided data for pathway-based analysis. \code{df_exp}: column 
# names -- sample id; row names -- gene symbols. Write exluded data to 
# "user_mutated_gene_exclude_pathway.txt" and 
# "user_expressed_gene_exclude_pathway.txt".
# 
# @param df_mut user-provided mutation data.
# @param df_exp user-provided expression data.
# @param sid_pattern A string regex.
# @return A list, first element is mutation data, second element is expression
#   data, third element is case/control sample index in expression data.
# @examples 
# 
# \dontrun{
# res <- preprocess_user_data_pathway(hnsc_mut, hnsc_exp, c("-01$", "-11$"))
# }
# @export
preprocess_user_data_pathway <- function(df_mut, df_exp, sid_pattern) {
  ret <- vector(mode = "list", length = 3)
  # check user's data column name
  mut_colnames <- colnames(df_mut)
  if (length(mut_colnames) != 7|
      mut_colnames[1] != "Hugo_Symbol"|
      mut_colnames[2] != "Chromosome"|
      mut_colnames[3] != "Start_position"|
      mut_colnames[4] != "End_position"|
      mut_colnames[5] != "Reference_Allele"|
      mut_colnames[6] != "Tumor_Seq_Allele2"|
      mut_colnames[7] != "Tumor_Sample_Barcode") {
    stop("Your gene mutation data is invalid, please check again!", call. = FALSE)
  }
  
  mutated_symbols <- df_mut[[1]]
  expressed_symbols <- rownames(df_exp)
  mut_idx <- mutated_symbols %in% valid_kegg_symbols()
  exp_idx1 <- expressed_symbols %in% valid_TG_symbols(pos = TRUE)
  exp_idx2 <- expressed_symbols %in% valid_TG_symbols(pos = FALSE)
  exp_idx3 <- expressed_symbols %in% valid_TF_symbols(pos = TRUE)
  exp_idx4 <- expressed_symbols %in% valid_TF_symbols(pos = FALSE)
  exp_idx5 <- expressed_symbols %in% valid_kegg_symbols()
  exp_idx <- exp_idx1 | exp_idx2 | exp_idx3 | exp_idx4 | exp_idx5
  
  # a few genes are not analyzed and are logged to files.
  write.table(df_mut[!mut_idx, ], file = "user_mutated_gene_exclude_pathway.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(df_exp[!exp_idx, ], file = "user_expressed_gene_exclude_pathway.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  # drop some user's data
  mut_included <- df_mut[mut_idx, ]
  exp_included <- df_exp[exp_idx, ]
  
  case_pattern <- sid_pattern[1]
  control_pattern <- sid_pattern[2]
  
  mut_sid <- unique(mut_included[[7]])
  mut_case <- grepl(case_pattern, mut_sid)
  if (!all(mut_case)) {
    stop("Your mutation data contains invalid sample IDs, please check again!", call. = FALSE)
  }
  
  exp_sid <- colnames(exp_included)
  exp_case <- grepl(case_pattern, exp_sid)
  exp_control <- grepl(control_pattern, exp_sid)
  if (!any(exp_case) | !any(exp_control)) {
    stop("Your expression data does not contain case/control data, please check again!", call. = FALSE)
  }
  
  if (sum(exp_case | exp_control) != length(exp_sid)) {
    stop("Your expression data contains invalid sample IDs, please check again!", call. = FALSE)
  }
  
  if (!setequal(mut_sid, exp_sid[exp_case])) {
    warning("Your mutation data's sample IDs probably does not match your expression data's sample IDs, please take carefull note of the analysis result!", call. = FALSE)
  }
  
  exp_case_sid_sorted <- sort(exp_sid[exp_case])
  exp_control_sid_sorted <- sort(exp_sid[exp_control])
  ret_exp_case_idx <- match(exp_case_sid_sorted, exp_sid)
  ret_exp_control_idx <- match(exp_control_sid_sorted, exp_sid)
  
  # return
  ret[[1]] <- mut_included
  ret[[2]] <- exp_included
  ret[[3]] <- list(case = ret_exp_case_idx, control = ret_exp_control_idx)
  ret
}

# Get valid transcription factor gene symbols.
# 
# Get valid transcription factor gene symbols.
# 
# @param pos Logical value.
# @return A character vector of valid transcription factor gene symbols.
# @examples 
# 
# \dontrun{
# valid_TF_symbols()
# }
# @export
valid_TF_symbols <- function(pos = TRUE) {
  ret <- NULL
  if (pos) {
    ret <- TF_target_pos$TF
    ret <- unique(ret)
  } else {
    ret <- TF_target_neg$TF
    ret <- unique(ret)
  }
  ret
}
# Get valid target gene symbols.
# 
# Get valid target gene symbols.
# 
# @param pos Logical value.
# @return A character vector of valid target gene symbols.
# @examples 
# 
# \dontrun{
# valid_TG_symbols()
# }
# @export
valid_TG_symbols <- function(pos = TRUE) {
  ret <- NULL
  if (pos) {
    ret <- TF_target_pos$TG
    ret <- unique(ret)
  } else {
    ret <- TF_target_neg$TG
    ret <- unique(ret)
  }
  ret
}
# Get valid KEGG gene symbols.
# 
# Get valid KEGG gene symbols.
# 
# @return A character vector of valid KEGG gene symbols.
# @examples 
# 
# \dontrun{
# valid_kegg_symbols()
# }
# @export
valid_kegg_symbols <- function() {
  ret <- lapply(kegg_pathway, function(x) {
    names(igraph::V(x))
  })
  
  ret <- base::Reduce(base::union, ret)
  unique(ret)
}
# Get valid gene symbols in the PPI.
# 
# Get valid gene symbols in the PPI.
# 
# @return A character vector of valid transcription factor gene symbols.
# @examples 
# 
# \dontrun{
# valid_ppi_symbols()
# }
# @export
valid_ppi_symbols <- function() {
  ret <- names(igraph::V(net_template))
  unique(ret)
}

#' Example data
#' 
#' Example data
"hnsc_mut_part"
#' Example data
#' 
#' Example data
"hnsc_exp_part"
#' Example data
#' 
#' Example data
"hnsc_expressed_genes"