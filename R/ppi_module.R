# Greedy algorithm for maximum coverage problem.
# 
# Perform greedy algorithm for maximum coverage problem.
# 
# @param s A list for the problem.
# @return An optimization result of the problem.
# @examples 
# s <- setNames(list(
# c(1, 2, 3, 8, 9, 10),
# c(1, 2, 3, 4, 5),
# c(4, 5, 7),
# c(5, 6, 7),
# c(6, 7, 8, 9, 10)), 
# c("s1", "s2", "s3", "s4", "s5"))
# res <- greedy_set_cover(s)
# exp_res <- list(
# s1 = c(1, 2, 3, 8, 9, 10),
# s3 = c(4, 5, 7),
# s4 = c(6))
# @export
greedy_set_cover <- function(s) {
  stopifnot(is.list(s))
  
  s <- compact(s)
  ret <- list()
  u <- Reduce(base::union, s)
  
  s_copy <- s
  while(length(u) > 0) {

    l <- vapply(s_copy, length, FUN.VALUE = integer(1))
    # which.max return first max element index
    i <- which.max(l)
    
    # update c_idx, u, ss
    ret <- c(ret, s_copy[i])
    tmp <- s_copy[[i]]
    u <- setdiff(u, tmp)
    s_copy <- lapply(s_copy, function(x) setdiff(x, tmp))
    s_copy <- compact(s_copy)
  }
  
  ret
}

# Find optimal solution for the mutant genes and differentially expressed genes.
# 
# Find optimal solution for the mutant genes and differentially expressed genes.
# 
# @param x1 A character vector of mutant genes.
# @param  x2 A character vector of differentially expressed genes.
# @param nei_order The order of neighborhood from the mutant gene.
# @param method Algorithm used to solve the problem.
# @return A list for the solution.
# @export
find_sig_genes <- function(x1, x2, nei_order, method = greedy_set_cover) {
  net_gene <- igraph::V(net_template)$name
  v <- intersect(x1, net_gene)
  
  nei <- igraph::ego(net_template, nei_order, v)
  nei_v <- lapply(nei, function(x) intersect(x2, x$name))
  
  s <- setNames(nei_v, v)
  method(s)
}