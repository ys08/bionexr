#' Pipeline of PPI-based analysis.
#' 
#' This is the PPI-based analysis pipeline.
#' 
#' \code{df_user_mut} is a dataframe with 7 columns, which are "Hugo_Symbol", 
#' "Chromosome", "Start_position", "End_position", "Reference_Allele", 
#' "Tumor_Seq_Allele2", "Tumor_Sample_Barcode". \code{df_user_exp} is a 
#' dataframe, row names are gene symbols, and column names are sample IDs.
#' \code{sid_pattern} is a regex string that can distingish case ID and control
#' ID.
#' @param df_user_mut Mutation dataframe.
#' @param df_user_exp Gene expression dataframe.
#' @param sid_pattern A regex string to match case and control IDs.
#' @param mcutoff A float number between 0 and 6, default is 2.
#' @param fraction A float number between 0 and 1, default is 0.6.
#' @param nei_order The order of neighborhood from the mutant gene, default is 2.
#' @param mweight weight of mutation score for combined score, default is 0.6.
#' @param padj Adjusted p value threshold for differential analysis, default is 0.01.
#' @param log2_fold_change log2 transformed fold change threshold for 
#'   differential analysis, default is 2.
#' @param jobname unique job name.
#' @param raw_count Logical value, set TRUE to process RSEM-based expression 
#'   data.
#' @param use_cache  Logical value, set TRUE to save intermediate result.
#' @return A list.
#' @examples 
#' 
#' \dontrun{
#' mut_data <- firehose_get("HNSC", "mutation", run_date = "2015_08_21", run_type = "stddata")
#' mut_data <- mut_data[[1]]
#' mut_sample_ids <- unique(mut_data[[7]])
#' 
#' exp_data <- firehose_get("HNSC", "expression", run_date = "2015_08_21", run_type = "stddata")
#' exp_data <- exp_data[[1]]
#' exp_sample_ids <- colnames(exp_data)
#' 
#' common_case <- intersect(mut_sample_ids, exp_sample_ids)
#' exp_control <- grepl("-11$", exp_sample_ids)
#' hnsc_mut <- mut_data[mut_data[[7]] %in% common_case, ]
#' hnsc_exp <- exp_data[, (exp_sample_ids %in% common_case) | exp_control]
#' 
#' prepare_ma()
#' res <- perform_main_ppi(hnsc_mut, hnsc_exp, jobname = "HNSC", use_cache = TRUE)
#' }
#' @export
perform_main_ppi <- function(df_user_mut, df_user_exp, sid_pattern = c("-01$", "-11$"), mcutoff = 2, fraction = 0.6, nei_order = 2, mweight = 0.6, padj = 0.01, log2_fold_change = 2, jobname = "tcga", raw_count = TRUE, use_cache = FALSE) {
  intermediate_cache_filenames <- paste0(jobname, c("_df_mut_score_long_ppi.rds", "_df_mut_score_wide_ppi.rds", "_whole_fc_genes_ppi.rds"))
  is_intermediate_file_exist <- sapply(intermediate_cache_filenames, file.exists)
  
  gene.res <- vector(mode = "list", length = 3)
  if (!use_cache) {
    gene.res <- perform_gene_ppi(df_user_mut, df_user_exp, sid_pattern = sid_pattern, raw_count = raw_count)
  } else {
    if (all(is_intermediate_file_exist)) {
      gene.res[[1]] <- readRDS(intermediate_cache_filenames[1])
      gene.res[[2]] <- readRDS(intermediate_cache_filenames[2])
      gene.res[[3]] <- readRDS(intermediate_cache_filenames[3])
    } else {
      gene.res <- perform_gene_ppi(df_user_mut, df_user_exp, sid_pattern = sid_pattern, raw_count = raw_count)
      saveRDS(gene.res[[1]], intermediate_cache_filenames[1])
      saveRDS(gene.res[[2]], intermediate_cache_filenames[2])
      saveRDS(gene.res[[3]], intermediate_cache_filenames[3])
    }
  }
  
  ret <- perform_network_ppi(gene.res[[2]], gene.res[[3]], mcutoff = mcutoff, fraction = fraction, nei_order = nei_order, mweight = mweight, padj = padj, log2_fold_change = log2_fold_change)
  invisible(ret)
}

#' Perform gene module of PPI-based approach.
#' 
#' Perform gene module of PPI-based approach.
#' 
#' @param df_user_mut Mutation dataframe.
#' @param df_user_exp Gene expression dataframe.
#' @param sid_pattern A regex string to match case and control IDs.
#' @param raw_count Logical value, set TRUE to process RSEM-based expression 
#'   data.
#' @return A list.
#' @examples 
#' 
#' \dontrun{
#' 
#' prepare_ma()
#' res.gene <- perform_gene_ppi(hnsc_mut_part, hnsc_exp_part)
#' }
#' @export
perform_gene_ppi <- function(df_user_mut, df_user_exp, sid_pattern = c("-01$", "-11$"), raw_count = TRUE) {

  # preprocess user's data: check, drop
  df_clean <- preprocess_user_data_ppi(df_user_mut, df_user_exp, sid_pattern)
  df_mut <- df_clean[[1]]
  df_exp <- df_clean[[2]]
  idx <- df_clean[[3]]
  
  if (nrow(df_mut) == 0 | nrow(df_exp) == 0) {
    stop("Your data is invalid, please check again!", call. = FALSE)
  }
  whole_fc_genes <- de_analysis(df_exp[, idx[[1]]], df_exp[, idx[[2]]], raw_count = raw_count)
  df_mut_score_long <- score_mutation(df_mut)
  df_mut_score_wide <- reshape_long_to_wide(df_mut_score_long)
  
  invisible(list(df_mut_score_long, df_mut_score_wide, whole_fc_genes))
}

#' Perform network module of PPI-based approach.
#' 
#' Perform network module of PPI-based approach.
#' 
#' @param df_mut_score_wide A dataframe.
#' @param whole_fc_genes A list.
#' @param mcutoff A float number between 0 and 6, default is 2.
#' @param fraction A float number between 0 and 1, default is 0.6.
#' @param nei_order The order of neighborhood from the mutant gene, default is 2.
#' @param mweight weight of mutation score for combined score, default is 0.6.
#' @param padj Adjusted p value threshold for differential analysis, default is 0.01.
#' @param log2_fold_change log2 transformed fold change threshold for 
#'   differential analysis, default is 2.
#' @return A list.
#' @examples 
#' 
#' \dontrun{
#' 
#' prepare_ma()
#' res.gene <- perform_gene_ppi(hnsc_mut_part, hnsc_exp_part)
#' res.network <- perform_network_ppi(res.gene[[2]], res.gene[[3]])
#' }
#' @export
perform_network_ppi <- function(df_mut_score_wide, whole_fc_genes, mcutoff = 2, fraction = 0.6, nei_order = 2, mweight = 0.6, padj = 0.01, log2_fold_change = 2) {
  # filter
  padj_value <- whole_fc_genes$padj_value
  padj_value_idx <- padj_value <= padj
  padj_value_idx <- na_to_value(padj_value_idx, FALSE)
  log2_fc <- whole_fc_genes$log2_fc
  log2_fc_idx <- abs(log2_fc) >= log2_fold_change
  log2_fc_idx <- na_to_value(log2_fc_idx, FALSE)
  
  filter_idx <- padj_value_idx & log2_fc_idx
  whole_fc_genes <- whole_fc_genes[filter_idx, ]
  
  # core steps
  fc_genes <- whole_fc_genes$gene_name
  mut_genes <- choose_mut_genes(df_mut_score_wide, mcutoff, fraction)
  # greedy algrithms
  net_res <- ppi_analysis(mut_genes, fc_genes, nei_order)
  # a naive scoring foluma
  net_score <- unlist(lapply(net_res, length))
  # aggregate gene mutation score by sum each sample's score
  mut_score <- apply(df_mut_score_wide[mut_genes, , drop = FALSE], 1, function(x) {
    sum(x, na.rm = TRUE)
  })
  #combine two score
  score <- data.frame(mutation = mut_score, network = 1, stringsAsFactors = FALSE)
  score[names(net_score), "network"] <- net_score
  # normalization
  score[] <- log_normalization(score)
  score$combine <- mweight * score$mutation + (1 - mweight) * score$network
  score <- score[order(score$combine, decreasing = TRUE), ]
  
  score <- round(score, 2)
  # add fold change to result
  res.net <- vector(mode = "list", length = length(net_res))
  for (i in seq_len(length(net_res))) {
    tmp_genes <- net_res[[i]]
    tmp_fc <- whole_fc_genes[fc_genes %in% tmp_genes, "log2_fc"]
    tmp_df <- data.frame(gene_name = tmp_genes, log2_fc = tmp_fc, stringsAsFactors = FALSE)
    res.net[[i]] <- tmp_df
  }
  res.net <- setNames(res.net, names(net_res))
  invisible(list(score, res.net))
}

# Core functions used for PPI-based analysis.
# 
# Core functions used for PPI-based analysis.
#
# @param mut_genes A character vector of mutated gene symbols.
# @param fc_genes A character vector of fold change gene symbols.
# @param nei_order The order of neighborhood from the mutant gene.
# @return A list for the solution.
ppi_analysis <- function(mut_genes, fc_genes, nei_order) {
  N_mut <- length(mut_genes)
  N_fc <- length(fc_genes)
  message("Explaining ", N_mut, " mutated genes and ", N_fc, " differentially expressed genes using protein-protein network!")
  ret <- find_sig_genes(mut_genes, fc_genes, nei_order)
  invisible(ret)
}



