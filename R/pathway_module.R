# Find subgraph between mutated gene and transcription factor in the pathway.
# 
# Find subgraph between mutated gene and transcription factor in the pathway.
# 
# @param mutated_gene A string of gene symbol.
# @param tf A string of gene symbol.
# @param pathway A list of igraph object.
# @return A list of igraph object.
# @export
find_mutation_downstream_subnet <- function(mutated_gene, tf, pathway) {
  ret <- NULL
  if (igraph::are.connected(pathway, mutated_gene, tf)) {
    a <- igraph::subcomponent(pathway, mutated_gene, mode = "out")
    b <- igraph::subcomponent(pathway, tf, mode = "in")
    c <- intersect(names(a), names(b))
    
    ret <- igraph::induced_subgraph(pathway, c)
  }
  ret
}
# find_mutation_downstream_subnet <- function(mutated_gene, tf, pathway) {
#   #stop_watch <- Sys.time()
#   ret <- NULL
#   if (igraph::are.connected(pathway, mutated_gene, tf)) {
#     all_path <- igraph::all_simple_paths(pathway, from = mutated_gene, to = tf)
#     all_node <- lapply(all_path, names)
#     all_node <- Reduce(base::union, all_node)
#     ret <- igraph::induced_subgraph(pathway, all_node)
#   }
#   
#   #elapsed_time <- Sys.time() - stop_watch
#   #message(elapsed_time)
#   ret
# }

# Find shortest path between mutated gene and transcription factor in the pathway.
# 
# Find shortest path between mutated gene and transcription factor in the pathway.
# 
# @param mutated_gene A string of gene symbol.
# @param tf A string of gene symbol.
# @param pathway A list of igraph object.
find_mutation_downstream_shortest <- function(mutated_gene, tf, pathway) {
  branches <- igraph::shortest_paths(pathway, from = mutated_gene, to = tf)
  branches
}

# Fisher's exact test.
# 
# Perform fisher's exact test.
subnet_enrichment <- function(subnet_genes, up_genes, all_genes) {
  # cat(length(subnet_genes), " subnet genes, ", length(up_genes), " dyregulated genes, ", length(intersect(subnet_genes, up_genes)), " overlap genes\n")
  fisher.res <- fisher_test(gene_list = subnet_genes, pathway = up_genes, all_list = all_genes)
  fisher.res
}

# Fisher's exact test.
# 
# Perform fisher's exact test.
# 
# @export
fisher_test <- function(gene_list, pathway, all_list) {
  a <- length(intersect(gene_list, pathway))
  b <- length(gene_list) - a
  c <- length(pathway) - a
  d <- length(all_list) - a - b - c
  matrix <- matrix(c(a, c, b, d), nrow = 2)
  fisher.test(matrix, alternative = "greater")$p.value
}

# Identify active transcription factor by fisher's exact test.
# 
# Transcription factor associated with significant number of up-regulated targets is considered to be active.
# 
# @param tf_list A character vector of transcription factors.
# @param up_gene_list A character vector of up-regulated genes.
# @param all_gene_list A character vector of all genes.
# @return A list, first element is p value, second element is target genes.
# @export
find_active_tf <- function(tf_list, up_gene_list, all_gene_list) {
  all_gene_list <- base::union(all_gene_list, valid_TF_symbols())
  all_gene_list <- base::union(all_gene_list, valid_TG_symbols())
  all_gene_list <- base::union(all_gene_list, valid_kegg_symbols())
  #cat(length(all_gene_list))
  
  ret <- vapply(tf_list, function(tf) {
    tg_idx <- TF_target_pos$TF == tf
    tg <- TF_target_pos$TG[tg_idx]
    fisher_test(tg, up_gene_list, all_gene_list)
  }, FUN.VALUE = numeric(1))
  
  ret2 <- lapply(tf_list, function(tf) {
    tg_idx <- TF_target_pos$TF == tf
    tg <- TF_target_pos$TG[tg_idx]
    base::intersect(tg, up_gene_list)
  })
  names(ret2) <- tf_list
  # ret is fisher p value, ret2 is related target genes
  list(ret, ret2)
}

# tf associated with significant number of downregulated targets is considered to be inactive
find_inactive_tf <- function(tf_list, down_gene_list, all_gene_list) {
  ret <- vapply(tf_list, function(tf) {
    tg_idx <- TF_target_neg$TF == tf
    tg <- TF_target_neg$TG[tg_idx]
    fisher_test(tg, down_gene_list, all_gene_list)
  }, FUN.VALUE = numeric(1))
  ret
}

# Extract significant pathway branches from the pathway analysis result.
# 
# Extract significant pathway branches from the pathway analysis result.
# 
# @param pathway_result pathway analysis result.
# @param p_thres P value threshold
# @return A list of the pathway branches, p value <= \code{p_thres}
# @export
significant_branches_from_fisher_test <- function(pathway_result, p_thres) {
  ret <- NULL
  branch <- pathway_result[[1]]
  sig_matrix <- pathway_result[[2]]
  ret <- lapply(seq_along(branch), function(i) {
    p <- branch[[i]]
    m <- sig_matrix[[i]]
    
    if(length(p) == 0) {
      return(NULL)
    }
    
    for (mut_gene in rownames(m)) {
      for (tf_gene in colnames(m)) {
        if (m[mut_gene, tf_gene] >= p_thres) {
          p[[mut_gene]][[tf_gene]] <- NULL
        }
      }
      if (length(p[[mut_gene]]) == 0) {
        p[[mut_gene]] <- NULL
      }
    }
    
    if (length(p) == 0) {
      return(NULL)
    } else {
      return(p)
    }
  })
  
  names(ret) <- names(branch)
  ret[!unlist(lapply(ret, is.null))]
}


# library(igraph)

# setup graph
# g= graph.formula(A -+ B,
#                  A -+ C,
#                  B -+ C,
#                  B -+ D,
#                  B -+ E
# )
# plot(g, layout = layout.reingold.tilford(g, root="A"))