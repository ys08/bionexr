#' visualize PPI-based network.
#' 
#' Visualize PPI-based network, this function takes
#' \code{network_from_ppi} output as input.
#' 
#' @param net An igraph object.
#' @param ... Other parameters.
#' @examples 
#' 
#' \dontrun{
#' 
#' prepare_ma()
#' res.gene <- perform_gene_ppi(hnsc_mut_part, hnsc_exp_part)
#' res.network <- perform_network_ppi(res.gene[[2]], res.gene[[3]])
#' g <- network_from_ppi(res.network)
#' plot_ppi(g)
#' }
#' @export
plot_ppi <- function(net, ...) {
  igraph::plot.igraph(net, ...)
}

#' Generate network from PPI-based result.
#' 
#' Generate network from PPI-based result.
#' 
#' @param ppi_result PPI-based analysis result.
#' @param topn Default value is 10.
#' @return A igraph object.
#' @examples 
#' 
#' \dontrun{
#' 
#' prepare_ma()
#' res.gene <- perform_gene_ppi(hnsc_mut_part, hnsc_exp_part)
#' res.network <- perform_network_ppi(res.gene[[2]], res.gene[[3]])
#' g <- network_from_ppi(res.network)
#' }
#' @export
network_from_ppi <- function(ppi_result, topn = 10) {
  net.from <- character()
  net.to <- character()
  net.igraph <- NULL
  mut_genes <- rownames(ppi_result[[1]])[1:topn]
  up_genes <- character()
  down_genes <- character()
  
  for (mg in mut_genes) {
    mg_de_genes <- ppi_result[[2]][[mg]]
    if (is.null(mg_de_genes)) {
      net.from <- c(net.from, mg)
      net.to <- c(net.to, mg)
    } else {
      # find mutate gene's neighbour de genes
      up_genes <- c(up_genes, mg_de_genes$gene_name[mg_de_genes$log2_fc > 0])
      down_genes <- c(down_genes, mg_de_genes$gene_name[mg_de_genes$log2_fc < 0])
      for (de_g in mg_de_genes$gene_name) {
        net.from <- c(net.from, mg)
        net.to <- c(net.to, de_g)
      }
    }
  }
  
  net.df <- data.frame(from = net.from, to = net.to, stringsAsFactors = FALSE)
  net.df <- unique(net.df)
  node_names <- unique(c(net.df[, 1], net.df[, 2]))
  mut_idx <- node_names %in% mut_genes
  up_idx <- node_names %in% up_genes
  down_idx <- node_names %in% down_genes
  
  # node attribute dataframe for plot purpose
  n_node <- length(node_names)
  nodes <- data.frame(node_names, 
                      color = character(n_node), 
                      frame.color = character(n_node), 
                      label = character(n_node), 
                      label.cex = numeric(n_node), 
                      label.font = numeric(n_node), 
                      size = numeric(n_node), 
                      stringsAsFactors = FALSE)
  # node color
  nodes$color[up_idx] <- "red"
  nodes$color[down_idx] <- "blue"
  nodes$color[mut_idx] <- adjustcolor("white", 0)
  # node frame color
  nodes$frame.color = NA
  # node label
  nodes$label[mut_idx] <- node_names[mut_idx]
  nodes$label[!mut_idx] <- NA
  # node label size
  nodes$label.cex = 0.6
  # node label style, 2 for bold style
  nodes$label.font = 2
  # node size
  nodes$size <- 2
  
  net.igraph <- igraph::graph_from_data_frame(net.df, directed = FALSE, vertices = nodes)
  net.igraph <- igraph::simplify(net.igraph)
  invisible(net.igraph)
}


