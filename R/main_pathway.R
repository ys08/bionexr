#' Pipeline of pathway-based analysis.
#' 
#' This is the pathway-based analysis pipeline.
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
#' res <- perform_main_pathway(hnsc_mut, hnsc_exp, jobname = "HNSC", use_cache = TRUE)
#' }
#' @export
perform_main_pathway <- function(df_user_mut, df_user_exp, sid_pattern = c("-01$", "-11$"), mcutoff = 2, fraction = 0.6, padj = 0.01, log2_fold_change = 2, jobname = "tcga", raw_count = TRUE, use_cache = FALSE) {
  
  intermediate_cache_filenames <- paste0(jobname, c("_df_mut_score_long_pathway.rds", "_df_mut_score_wide_pathway.rds", "_whole_fc_genes_pathway.rds"))
  is_intermediate_file_exist <- sapply(intermediate_cache_filenames, file.exists)
  
  gene.res <- vector(mode = "list", length = 3)
  if (!use_cache) {
    gene.res <- perform_gene_pathway(df_user_mut, df_user_exp, sid_pattern = sid_pattern, raw_count = raw_count)
  } else {
    if (all(is_intermediate_file_exist)) {
      gene.res[[1]] <- readRDS(intermediate_cache_filenames[1])
      gene.res[[2]] <- readRDS(intermediate_cache_filenames[2])
      gene.res[[3]] <- readRDS(intermediate_cache_filenames[3])
    } else {
      gene.res <- perform_gene_pathway(df_user_mut, df_user_exp, sid_pattern = sid_pattern, raw_count = raw_count)
      saveRDS(gene.res[[1]], intermediate_cache_filenames[1])
      saveRDS(gene.res[[2]], intermediate_cache_filenames[2])
      saveRDS(gene.res[[3]], intermediate_cache_filenames[3])
    }
  }
  
  # identify expressed genes using user's orignal data
  expressed_genes <- identify_expressed_genes(df_user_exp, rsem_thres = 1, frac_thres = 0.5)
  
  ret <- perform_network_pathway(gene.res[[2]], gene.res[[3]], expressed_genes, mcutoff = mcutoff, fraction = fraction, padj = padj, log2_fold_change = log2_fold_change)
  invisible(ret)
  
}

#' Perform gene module of pathway-based approach.
#' 
#' Perform gene module of pathway-based approach.
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
#' res.gene <- perform_gene_pathway(hnsc_mut_part, hnsc_exp_part)
#' }
#' @export
perform_gene_pathway <- function(df_user_mut, df_user_exp, sid_pattern = c("-01$", "-11$"), raw_count = TRUE) {
  # preprocess user's data: check, drop
  df_clean <- preprocess_user_data_pathway(df_user_mut, df_user_exp, sid_pattern)
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

#' Perform network module of pathway-based approach.
#' 
#' Perform network module of pathway-based approach.
#' 
#' @param df_mut_score_wide A dataframe.
#' @param whole_fc_genes A list.
#' @param expressed_genes A character vector.
#' @param mcutoff A float number between 0 and 6, default is 2.
#' @param fraction A float number between 0 and 1, default is 0.6.
#' @param padj Adjusted p value threshold for differential analysis, default is 0.01.
#' @param log2_fold_change log2 transformed fold change threshold for 
#'   differential analysis, default is 2.
#' @return A list.
#' @examples 
#' 
#' \dontrun{
#' 
#' prepare_ma()
#' res.gene <- perform_gene_pathway(hnsc_mut_part, hnsc_exp_part)
#' res.network <- perform_network_pathway(res.gene[[2]], res.gene[[3]], hnsc_expressed_genes)
#' }
#' @export
perform_network_pathway <- function(df_mut_score_wide, whole_fc_genes, expressed_genes, mcutoff = 2, fraction = 0.6, padj = 0.01, log2_fold_change = 2) {
  
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
  mut_genes <- choose_mut_genes(df_mut_score_wide, mcutoff, fraction)

  up_idx <- whole_fc_genes$log2_fc > 0
  up_genes <- whole_fc_genes$gene_name[up_idx]
  down_genes <- whole_fc_genes$gene_name[!up_idx]
  
  active_branch <- identify_active_branch(mut_genes, up_genes, whole_fc_genes$gene_name, expressed_genes)
  active_branch
  
}

# Core functions used for identify active pathway branches.
# 
# Core functions used for identify active pathway branches.
# 
# @param mut_genes A character vector of mutated gene symbols.
# @param up_genes A character vector of up-regulated gene symbols.
# @param fc_genes A character vector of fold change gene symbols.
# @param expressed_genes A character vector of expressed gene symbols.
# @return A list of active branches.
identify_active_branch <- function(mut_genes, up_genes, fc_genes, expressed_genes) {
  N_mut <- length(mut_genes)
  N_fc <- length(up_genes)
  #kegg_up_genes <- base::intersect(valid_kegg_symbols(), up_genes)
  #kegg_fc_genes <- base::intersect(valid_kegg_symbols(), fc_genes)
  all_gene_list <- base::union(expressed_genes, valid_kegg_symbols())
  # message("Total ", length(all_gene_list), " genes, explaining ", N_mut, " mutated genes and ", N_fc, " up-regulated genes using KEGG pathway!")
  #message(length(all_gene_list), " KEGG genes.")
  #message(length(kegg_fc_genes), " fold-change genes are in KEGG genes.")
  #message(length(kegg_up_genes), " up-regulated genes are in KEGG genes.")

  expressed_tf <- base::intersect(valid_TF_symbols(), expressed_genes)
  
  active_tf_res <- find_active_tf(expressed_tf, up_genes, expressed_genes)
  active_tf <- active_tf_res[[1]]
  active_idx <- (active_tf <= 0.05)
  # cat(sum(active_idx))
  active_tf <- expressed_tf[active_idx]
  active_tf_target <- active_tf_res[[2]][active_idx]
  #cat(length(active_tf_target))
  #cat(sum(active_tf<=0.05))
  # kegg_pathway is internal data
  internal_pathway <- kegg_pathway
  # find pathway that enrich up-regulated genes
  fisher_pathway <- lapply(internal_pathway, function(p) {
    fisher_test(names(igraph::V(p)), up_genes, all_gene_list)
  })
  enrich_pathway <- internal_pathway[fisher_pathway <= 0.05]
  # ===
  ret.branch <- vector(mode = "list", length = length(internal_pathway))
  ret.significance <- vector(mode = "list", length = length(internal_pathway))
  names(ret.branch) <- names(internal_pathway)
  names(ret.significance) <- names(internal_pathway)
  
  f_log <- "identify_active_branch.log.txt"
  write(as.character(Sys.time()), file = f_log)
  write(paste0(paste0(names(enrich_pathway), collapse = ", "), " enrich up-regulated genes."), file = f_log, append = TRUE)
  for (i in seq_along(internal_pathway)) {
    p_i <- internal_pathway[[i]]
    p_node_symbols <- names(igraph::V(p_i))
    p_mut <- base::intersect(mut_genes, p_node_symbols)
    p_tf <- base::intersect(active_tf, p_node_symbols)
    
    branch <- vector(mode = "list", length = length(p_mut))
    names(branch) <- p_mut
    write(names(internal_pathway)[i], file = f_log, append = TRUE)
    write(paste0(length(p_mut), " * ", length(p_tf)), file = f_log, append = TRUE)
    write(paste0("MUT GENES: ", paste0(p_mut, collapse = "\t")), file = f_log, append = TRUE)
    write(paste0("ACTIVE TF: ", paste0(p_tf, collapse = "\t")), file = f_log, append = TRUE)
    if (length(p_mut) == 0 | length(p_tf) == 0) {
      branch <- list()
    } else {
      for (g_mut in p_mut) {
        tmp <- vector(mode = "list", length = length(p_tf))
        names(tmp) <- p_tf
        for (g_tf in p_tf) {
          # write(paste0(g_mut, " ", g_tf), file = f_log, append = TRUE)
          tmp[[g_tf]] <- find_mutation_downstream_subnet(g_mut, g_tf, p_i)
        }
        branch[[g_mut]] <- tmp
      }
    }
    # significant test for each branch
    sig_matrix <- matrix(nrow = length(p_mut), ncol = length(p_tf))
    rownames(sig_matrix) <- p_mut
    colnames(sig_matrix) <- p_tf
    
    if (length(p_mut) == 0 | length(p_tf) == 0) {
      sig_matrix <- list()
    } else {
      for (g_mut in p_mut) {
        for(g_tf in p_tf) {
          p_b <- branch[[g_mut]][[g_tf]]
          if (igraph::is_igraph(p_b)) {
            subnet_genes <- names(igraph::V(p_b))
            subnet_genes <- setdiff(subnet_genes, g_mut)
            fisher_res <- subnet_enrichment(subnet_genes, up_genes, all_gene_list)
            sig_matrix[g_mut, g_tf] <- fisher_res
          } else {
            sig_matrix[g_mut, g_tf] <- Inf
          }
        }
      }
    }

    ret.branch[[i]] <- branch
    ret.significance[[i]] <- sig_matrix
    write(as.character(Sys.time()), file = f_log, append = TRUE)
  }
  write(as.character(Sys.time()), file = f_log, append = TRUE)
  list(branch = ret.branch, significance = ret.significance, tf_target = active_tf_target)
}



