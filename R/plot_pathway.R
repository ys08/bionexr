#' visualize pathway-based network.
#' 
#' Visualized pathway-based network, this function takes
#' \code{network_from_significant_branches} output as input.
#' 
#' @param net An igraph object.
#' @param ... Other parameters.
#' @examples 
#' 
#' \dontrun{
#' 
#' prepare_ma()
#' res.gene <- perform_gene_pathway(hnsc_mut_part, hnsc_exp_part)
#' res.network <- perform_network_pathway(res.gene[[2]], res.gene[[3]], hnsc_expressed_genes)
#' 
#' g <- network_from_significant_branches(res.network)
#' plot_pathway(g)
#' }
#' @export
plot_pathway <- function(net, ...) {
  res <- layout_by_vertex_attr(net)
  layout <- res[[1]]
  
  igraph::plot.igraph(net, layout = layout, edge.arrow.size = 0.5, ...)
  abline(a = res[[2]], b = 0)
  abline(a = res[[3]], b = 0)
}

#' Generate network from pathway branches.
#' 
#' Generate network from pathway branches.
#' 
#' @param pathway_result Pathway-based analysis result.
#' @param p_thres P value threshold, default 0.05.
#' @return A igraph object.
#' @examples 
#' 
#' \dontrun{
#' 
#' prepare_ma()
#' res.gene <- perform_gene_pathway(hnsc_mut_part, hnsc_exp_part)
#' res.network <- perform_network_pathway(res.gene[[2]], res.gene[[3]], hnsc_expressed_genes)
#' 
#' g <- network_from_significant_branches(res.network)
#' }
#' @export
network_from_significant_branches <- function(pathway_result, p_thres = 0.05) {
  sig_branches <- significant_branches_from_fisher_test(pathway_result, p_thres)
  tf_tg_table <- pathway_result[[3]]
  
  net.from <- character()
  net.to <- character()
  net.igraph <- NULL
  mut_genes <- character()
  tf_genes <- character()
  
  for ( i in seq_along(sig_branches)) {
    p <- sig_branches[[i]]
    p_mut_genes <- names(p)
    mut_genes <- c(mut_genes, p_mut_genes)
    
    for (p_mg in p_mut_genes) {
      mg_tfs <- names(p[[p_mg]])
      tf_genes <- c(tf_genes, mg_tfs)
      
      for (tf in mg_tfs) {
        net.from <- c(net.from, p_mg)
        net.to <- c(net.to, tf)
        
        tf_tgs <- tf_tg_table[[tf]]
        for (tg in tf_tgs) {
          net.from <- c(net.from, tf)
          net.to <- c(net.to, tg)
        }
      }
    }
  }
  
  mut_genes <- unique(mut_genes)
  tf_genes <- unique(tf_genes)
  
  net.df <- data.frame(from = net.from, to = net.to, stringsAsFactors = FALSE)
  net.df <- unique(net.df)
  
  node_names <- unique(c(net.df[, 1], net.df[, 2]))
  mut_type_idx <- node_names %in% mut_genes
  tf_type_idx <- node_names %in% tf_genes
  tg_types_idx <- !(mut_type_idx | tf_type_idx)
  
  # node attribute dataframe for plot purpose
  n_node <- length(node_names)
  nodes <- data.frame(node_names, 
                      level = numeric(n_node),
                      color = character(n_node),
                      frame.color = character(n_node), 
                      label = character(n_node), 
                      label.cex = numeric(n_node), 
                      label.font = numeric(n_node), 
                      label.dist = numeric(n_node), 
                      size = numeric(n_node), 
                      stringsAsFactors = FALSE)
  # set node level
  nodes$level[tg_types_idx] = 3
  nodes$level[tf_type_idx] = 2
  nodes$level[mut_type_idx] = 1
  # set node color
  nodes$color[tg_types_idx] = "green"
  nodes$color[tf_type_idx] = "blue"
  nodes$color[mut_type_idx] = "red"
  # node frame color
  nodes$frame.color = NA
  # node label
  nodes$label[mut_type_idx] <- node_names[mut_type_idx]
  nodes$label[tf_type_idx] <- node_names[tf_type_idx]
  nodes$label[tg_types_idx] <- NA
  # node label size
  nodes$label.cex = 0.6
  # node label distance
  nodes$label.dist = 0.6
  # node label style, 2 for bold style
  nodes$label.font = 2
  # node size
  nodes$size <- 5
  
  net.igraph <- igraph::graph_from_data_frame(net.df, directed = TRUE, vertices = nodes)
  
  invisible(net.igraph)
}

# Generate layout by vertex attribute.
# 
# Generate a layout that group vertex with the same attribute.
# 
# @param obj_igraph A igraph object.
# @param n Size of the grid, default is 1000.
# @param n_row1 Number of rows of level 1, default is 50.
# @param n_row2 Number of rows of level 2, default is 50.
# @param n_sep Number of blank between levels, default is 200.
# @return An igraph layout object.
# @export
layout_by_vertex_attr <- function(obj_igraph, n = 1000, n_row1 = 50, n_row2 = 50, n_sep = 200) {  
  node_levels <- igraph::V(obj_igraph)$level
  mut_nodes <- which(node_levels == 1)
  tf_nodes <- which(node_levels == 2)
  tg_nodes <- which(node_levels == 3)
  
  n_row3 <- n - n_row1 - n_row2
  mut_grid <- matrix(rep(0, n * n_row1), ncol = n, byrow = TRUE)
  tf_grid <- matrix(rep(0, n * n_row2), ncol = n, byrow = TRUE)
  tg_grid <- matrix(rep(0, n * n_row3), ncol = n, byrow = TRUE)
  
  mut_grid <- random_location_for_few(mut_grid, mut_nodes)
  tf_grid <- random_location_for_few(tf_grid, tf_nodes)
  tg_grid <- random_location_for_many(tg_grid, tg_nodes)
  sep_grid <- matrix(rep(0, n * n_sep), ncol = n)
  
  grid <- rbind(mut_grid, sep_grid, tf_grid, sep_grid, tg_grid)
  l <- grid_layout(grid)
  
  # calculate line position
  l1 <- 1 - 2 * (n_row1 + n_sep / 2) / nrow(grid)
  l2 <- 1 - 2 * (n_row1 + n_row2 + 1.5 * n_sep) / nrow(grid)
  # cat(l1, l2)
  return(list(l, l1, l2))
}

random_location_for_few <- function(m, v) {
  i <- floor(ncol(m) / length(v))
  loc_x <- sample(seq(from = 1, to = ncol(m), by = i), length(v))
  loc_y <- sample(nrow(m), length(v), replace = TRUE)
  
  for (j in seq_along(loc_x)) {
    m[loc_y[j], loc_x[j]] <- v[j]
  }
  m
}

random_location_for_many <- function(m, v) {
  loc <- sample(nrow(m) * ncol(m), length(v))
  m[loc] <- v
  m
}

grid_layout <- function(x) {
  LmatX <- seq(-1, 1, length = ncol(x))
  LmatY <- seq(1, -1, length = nrow(x))
  
  loc <- t(sapply(1: max(x), function(y) which(x == y, arr.ind = T)))
  layout <- cbind(LmatX[loc[, 2]], LmatY[loc[, 1]])
  return(layout)
}