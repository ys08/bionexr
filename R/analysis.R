# Perform concensus clustering.
# 
# Perform concensus clustering.
# 
# Perform concensus clustering, user can provide the number of clusters in
# runtime. \code{df_exp} row names are genes, and column names are sample IDs.
# @param df_exp A dataframe of gene expression data. 
# @param topn Number of genes for clustering.
# @return A dataframe: sample ID and sample class
# @examples 
# 
# \dontrun{
# res <- cc_analysis(hnsc_exp)
# }
# @export
# cc_analysis <- function(df_exp, topn = 2000) {
#   if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
#     stop("This function needs ConsensusClusterPlus libraries installed, albeit they cannot be found. Check out the installation instructions!\n", 
#          call. = FALSE)
#   }
#   # log transformation
#   df_exp <- log2(df_exp + 1)
#   # select topn variant genes
#   mad_value <- apply(df_exp, 1, mad)
#   if (topn > nrow(df_exp)) {
#     topn = nrow(df_exp)
#   }
#   df_exp_subset <- df_exp[order(mad_value, decreasing = T)[1:topn], ]
#   df_exp_subset <- sweep(df_exp_subset, 1, apply(df_exp_subset, 1, median, na.rm = T))
#   # perform consensus cluster
#   message("Consensus cluster parameters: maxK = 20, reps = 1000, pItem = 0.8, clusterAlg = \"hc\", distance = \"pearson\", plot = \"pdf\"")
#   cc_result <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(df_exp_subset), 
#                                                           maxK = 20, 
#                                                           reps = 1000, 
#                                                           pItem = 0.8, 
#                                                           clusterAlg = "hc", 
#                                                           distance = "pearson",
#                                                           plot = "pdf")
#   # user specify the best cluster
#   k <- read_integer("Enter the kth [1-20] cluster you selected : ")
#   named_vector <- cc_result[[k]]$consensusClass
#   ret <- data.frame(sample_id = names(named_vector), 
#                     sample_class = paste0("g", named_vector),
#                     stringsAsFactors = FALSE)
#   
#   row.names(ret) <- NULL
#   invisible(ret)
# }

# Perform differential expression analysis.
# 
# Perform differential expression analysis.
# 
# Perform differential expression analysis. For RSEM-based data, use "DESeq2"
# package to test significant differentially expressed genes. For protein-based
# data, use t test instead.
# 
# @param df_case_exp A dataframe of case sample expression data.
# @param df_control_exp A dataframe of control sample expression data.
# @param raw_count Set to TRUE if the data is RSEM based, FALSE if the data is 
#   protein based.
# @return A dataframe with 3 columns: gene name, log2 fold change value, 
#   adjusted p value.
# @examples 
# 
# \dontrun{
# res <- de_analysis(hnsc_exp[, case_idx], hnsc_exp[, control_idx], 0.05, 2)
# }
# @export
de_analysis <- function(df_case_exp, df_control_exp, raw_count = TRUE) {
  ret <- NULL
  
  if (raw_count) {
    ret <- de_analysis_rna(df_case_exp, df_control_exp)
  } else {
    ret <- de_analysis_protein(df_case_exp, df_control_exp)
  }
  
  invisible(ret)
}

# 20501*326 matrix takes 2537.92 seconds to run
de_analysis_rna <- function(df_case_exp, df_control_exp) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("This function needs DESeq2 library installed, albeit it cannot be found. Check out the installation instructions!\n", 
         call. = FALSE)
  }
  
  df_exp <- cbind(df_case_exp, df_control_exp)
  
  sample_ids <- colnames(df_exp)
  metadata <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  metadata$condition <- NA
  metadata[metadata$sample_id %in% colnames(df_case_exp), "condition"] <- "case"
  metadata[metadata$sample_id %in% colnames(df_control_exp), "condition"] <- "con"
  metadata$condition <- as.factor(metadata$condition)
  
  # diff analysis
  df_exp <- round(df_exp)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = df_exp, colData = metadata, design = ~condition)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  
  gene_name <- rownames(res)
  padj_value <- res$padj
  log2_fc <- res$log2FoldChange
  
  ret <- data.frame(gene_name = gene_name, 
                    log2_fc = log2_fc, 
                    padj_value = padj_value,
                    stringsAsFactors = FALSE)
  invisible(ret)
}

de_analysis_protein <- function(df_case_exp, df_control_exp) {
  log2_case_exp <- log2(df_case_exp + 1)
  log2_control_exp <- log2(df_control_exp + 1)
  
  res.p_value <- vector(mode = "numeric", length = nrow(log2_case_exp))
  res.log2_fold_change <- vector(mode = "numeric", length = nrow(log2_case_exp))
  for (i in seq_len(nrow(log2_case_exp))) {
    res <- t.test(log2_case_exp[i, ], log2_control_exp[i, ])
    res.p_value[i] <- res$p.value
    res.log2_fold_change[i] <- res$estimate[1] - res$estimate[2]
  }
  res.padj <- p.adjust(res.p_value, method = "BH")
  
  ret <- data.frame(gene_name = row.names(df_case_exp), 
                    log2_fc = res.log2_fold_change, 
                    padj_value = res.padj,
                    stringsAsFactors = FALSE)
  invisible(ret)
}

# fisher's exact test for enrichment of gene list
# fomula from http://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
# x: test list; g: gene list; topn: test topn genes in test list 
fisher_enrichment_test <- function(x, g, topn = 50) {
  r1 <- x[1:topn]
  r2 <- x[(topn + 1):length(x)]
  r11 <- length(intersect(r1, g))
  r21 <- length(intersect(r2, g))
  r12 <- length(r1) - r11
  r22 <- length(r2) - r21
  
  m <- matrix(c(r11, r12, r21, r22), nrow = 2, byrow = TRUE)
  res.fisher <- fisher.test(m, alternative = "greater")
  
  res.confusion <- confusion_test(x, g, topn)
  
  ret <- c(enrichment_p_value = round(res.fisher$p.value, 2), res.confusion[c("sensitivity", "specificity", "fdr")])
  ret
}

# evaluate result, return a confusion matrix
confusion_test <- function(x, g, topn = 50) {
  positive <- intersect(g, x)
  negtive <- setdiff(x, positive)
  
  predict_positive <- x[1:topn]
  predict_negtive <- x[(topn + 1):length(x)]
  
  tp <- length(intersect(positive, predict_positive))
  tn <- length(intersect(negtive, predict_negtive))
  fp <- length(intersect(negtive, predict_positive))
  fn <- length(intersect(positive, predict_negtive))
  
  n <- 2
  sensitivity <- tp / (tp + fn)
  sensitivity <- round(sensitivity, n)
  specificity <- tn / (tn + fp)
  specificity <- round(specificity, n)
  precision <- tp / (tp + fp)
  precision <- round(precision, n)
  accuracy <- (tp + tn) / (length(positive) + length(negtive))
  accuracy <- round(accuracy, n)
  fdr <- fp / (tp + fp)
  fdr <- round(fdr, n)
  
  ret <- c(sensitivity = sensitivity,
           specificity = specificity,
           precision = precision,
           accuracy = accuracy,
           fdr = fdr)
  
  ret
}

#' Identify expressed genes.
# 
#' Identify expressed genes based on user-provided expression data.
#' 
#' @param df_user_exp A dataframe.
#' @param rsem_thres Default value is 1.
#' @param frac_thres Default value is 0.5.
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
#' hnsc_exp <- exp_data[, (exp_sample_ids %in% common_case) | exp_control]
#' expressed_genes <- identify_expressed_genes(hnsc_exp, rsem_thres = 1, frac_thres = 0.5)
#' }
#' @export
identify_expressed_genes <- function(df_user_exp, rsem_thres = 1, frac_thres = 0.5) {
  ret <- NULL
  
  N <- ncol(df_user_exp)
  n_thres <- ceiling(frac_thres * N)
  df_logical <- df_user_exp >= rsem_thres
  n <- rowSums(df_user_exp)
  
  expressed_idx <- (n >= n_thres)
  
  ret <- rownames(df_user_exp)[expressed_idx]
  
  ret
}