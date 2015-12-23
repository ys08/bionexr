# Compute a function's execution time.
# 
# Wrap a function, record its execution time.
# 
# @param f Function to be wrapped.
# @return The return value of f.
# @examples 
# 
# \dontrun{
# timeit(f)(x)
# }
# @export
timeit <- function(f) {
  var_name <- deparse(substitute(f))
  function(...) {
    ptm <- proc.time()
    res <- f(...)
    elapsed_time <- (proc.time() - ptm)[[3]]
    elapsed_time <- round(elapsed_time, 2)
    message(var_name, " takes ", elapsed_time, " seconds to run\n")
    res
  }
}

# Defaults for NA/NULL values.
# 
# Defaults for NA/NULL values.
# 
# @param a One value.
# @param b Alternative value.
# @return The defaults for NA/NULL.
# @examples 
# 
# \dontrun{
# NA %||% 1
# NULL %||% 1
# }
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

# Remove NULLs from a list.
# 
# Remove NULLs from a list.
# 
# @param x A list.
# @return A list without NULL element.
# @examples 
# 
# \dontrun{
# l <- list(1, NULL, 2)
# compact(l)
# }
compact <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

# Print to console every n invocations to check on a long running process.
# 
# Print to console every n invocations to check on a long running process.
# 
# @param n An integer.
# @param f A function.
# @return The return value of f.
dot_every <- function(n, f) {
  i <- 1
  function(...) {
    if (i %% n == 0) cat(".")
    i <<- i + 1
    f(...)
  }
}

# Obtain local file size.
# 
# Obtain local file size.
# 
# @param file_name A string of local file name.
# @return The byte size of the file.
get_local_file_size <- function(file_name) {
  stopifnot(file.exists(file_name))
  file.size(file_name)
}

# Convert byte to human readable format.
# 
# Convert byte to human readable format.
#
# @param x A integer in byte.
# @return A string equivalent to x.
# @examples 
# 
# \dontrun{
# human_readable(10240)
# }
human_readable <- function(x) {
  a <- c(c(1, 1024, 1048576, 1073741824))
  b <- c('B', 'KB', 'MB', 'GB')
  
  idx <- sum(x >= a)
  paste0(round(x/a[idx], 1), b[idx])
}

# Default value for NA in a atomic vector.
# 
# Default value for NA in a atomic vector.
# 
# @param v An atomic vector.
# @param value Alternative value for NA.
# @return An atomic vector.
# @examples 
# 
# \dontrun{
# v <- c(1, NA, 2)
# na_to_value(v, 3)
# }
na_to_value <- function(v, value) {
  v[is.na(v)] <- value
  v
}

# Read integer from stdin.
# 
# Read integer from stdin.
# 
# @param msg A string.
read_integer <- function(msg = "Enter an integer: ") { 
  n <- readline(prompt = msg)
  if(!grepl("^[0-9]+$", n)) {
    return(read_integer(msg))
  }
  
  return(as.integer(n))
}

# Log normalize data.
# 
# Perform log normalization.
# 
# @param x A list.
# @return A normalized list.
# @examples 
# 
# \dontrun{
# l <- list(c1 = c(1, 2), c2 = c(3, 4))
# log_normalization(l)
# }
log_normalization <- function(x) {
  ret <- lapply(x, function(y) {
    M <- log10(max(y))
    log10(y)/M
  })
  
  ret
}

# Min-max normalization.
# 
# Perform min-max normalization.
# 
# @param x A list.
# @return A normalized list.
# @examples 
# 
# \dontrun{
# l <- list(c1 = c(1, 2), c2 = c(3, 4))
# min_max_normalization(l)
# }
min_max_normalization <- function(x) {
  ret <- lapply(x, function(y) {
    max_value <- max(y)
    min_value <- min(y)
    M <- max_value - min_value
    (y - min_value) / M
  })
  
  ret
}

# Z-score normalization.
# 
# Perform z-score normalization.
# 
# @param x A numeric matrix.
# @return A normalized numeric matrix.
# @examples 
# 
# \dontrun{
# m <- matrix(c(1, 2, 3, 4), nrow = 2)
# z_score_normalization(m)
# }
z_score_normalization <- function(x) {
  ret <- scale(x)
  ret
}

# Remove rows with too many NAs in a dataframe.
# 
# Remove rows with too many NAs in a dataframe.
# 
# @param X A dataframe.
# @param ratio A float number between 0 and 1.
# @return The row indexes to be removed.
# @examples 
# 
# \dontrun{
# d <- data.frame(replicate(5, c(NA, 1, NA, 3)))
# remove_low_coverage(d, 0.5)
# }
remove_low_coverage <- function(X, ratio) {
  count <- rowSums(is.na(X))
  id <- which(count >= ratio * ncol(X))
  
  return(id)
}

# Impute missing data.
# 
# Impute missing data.
# 
# @param X A dataframe.
# @param ratio A float number between 0 and 1.
# @param ... Additional parameters.
# @return A dataframe.
# impute_knn = function(X, ratio = 0.5, ...) {
#   if(!requireNamespace("impute", quietly = TRUE)) {
#     stop("This function needs the impute library installed, albeit it cannot be found. Check out the installation instructions!\n", 
#          call. = FALSE)
#   }
#   
#   id.col = remove_low_coverage(X, ratio)
#   id.row = remove_low_coverage(t(X), ratio)
#   
#   if (length(id.col) > 0) {
#     X = X[, -id.col]
#   }
#   if (length(id.row) > 0) {
#     X = X[-id.row, ]
#   }
#   
#   X = impute::impute.knn(data = X, ...)
#   X = t(X$data)
#   
#   return(X)
# }
