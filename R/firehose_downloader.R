#' Main function for downloading data from firehose.
#' 
#' Main function for downloading data from firehose.
#' 
#' @param disease_cohort use \code{list_cohorts} to check available disease cohorts.
#' @param data_type "mutation", "expression", or "Aggregate_AnalysisFeatures".
#' @param run_type use \code{list_runs} to check available run types.
#' @param run_date use \code{list_runs} to check valid run dates.
#' @param resume set TRUE to resume broken downloads.
#' @return R object for firehose data.
#' @examples 
#' 
#' \dontrun{
#' mut_data <- firehose_get("HNSC", "mutation", run_date = "2015_08_21", run_type = "stddata")
#' exp_data <- firehose_get("HNSC", "expression", run_date = "2015_08_21", run_type = "stddata")
#' }
#' @export
firehose_get <- function(disease_cohort, data_type , run_type = "stddata", run_date = "latest", resume = TRUE) {
  if (!all(requireNamespace("RCurl", quietly = TRUE), requireNamespace("XML", quietly = TRUE))) {
    stop("This function needs RCurl and XML libraries installed, albeit they cannot be found. Check out the installation instructions!\n", 
         call. = FALSE)
  }
  
  runs <- list_runs()
  cohorts <- list_cohorts()
  
  if (!all(disease_cohort %in% cohorts)) {
    err_msg <- "Incorrect disease cohorts, please run list_cohorts() to check available cohorts"
    stop(err_msg, call. = FALSE)
  }
  
  which_run <- paste(run_type, run_date, sep = "__")
  # map latest to exact date, such as xxxx_xx_xx
  if (run_date == "latest") {
    idx <- grepl(paste0("^", run_type), runs)
    
    if (!any(idx)) {
      err_msg <- "Incorrect runs, please run list_runs() to check available runs"
      stop(err_msg, call. = FALSE)
    }
    
    latest_run <- runs[idx][sum(idx)]
    message("Mapping ", which_run, " to ", latest_run)
    which_run <- latest_run
  }
  
  if (!any(which_run %in% runs)) {
    err_msg <- "Incorrect runs, please run list_runs() to check available runs"
    stop(err_msg, call. = FALSE)
  }
  
  # construct data url
  base_url <- "http://gdac.broadinstitute.org/runs"
  l1 <- which_run
  l2 <- "data"
  l3 <- disease_cohort
  l4 <- gsub("[a-z]|_", "", which_run)
  
  data_folder_url <- paste(base_url, l1, l2, l3, l4, sep = "/")
  # 301 moved permanently
  data_folder_url <- paste0(data_folder_url, "/")
  data_links <- list_data_file(data_folder_url)
  full_data_type <- paste(run_type, data_type, sep = "_")
  data_link <- match_data_type(data_links, full_data_type)
  
  # 1 data type per run
  dn_info <- lapply(data_link, function(x) download_file(x, resume = resume))
  downloaded_file <- unlist(lapply(dn_info, "[[", 3))
  
  # construct R data structure
  ret <- lapply(downloaded_file, function(x) dataframe_from_firehose(x, full_data_type))
  names(ret) <- downloaded_file
  invisible(ret)
}

# List available data files on firehose website.
# 
# List available data files on firehose website.
# 
# @param url firehose url.
# @return firehose data urls.
list_data_file <- function(url) {
  if (RCurl::url.exists(url)) {
    h <- RCurl::basicTextGatherer()
    res_code <- RCurl::curlPerform(url = url, writefunction = h$update)
    content <- h$value()
    # parse XML
    doc <- XML::htmlParse(content, asText = TRUE)
    
    file_name <- XML::xpathSApply(doc, "//table/..//a", XML::xmlGetAttr, "href")
    # real data with suffix tar.gz
    idx <- grepl("tar\\.gz$", file_name)
    
    data_link <- paste0(url, file_name[idx])
    data_link
  } else {
    stop(paste(url, "does not exist"), call. = FALSE)
  }
}

#' List available cohorts on firehose website.
#' 
#' List available cohorts on firehose website.
#' 
#' @return A character vector of cohorts.
#' @examples 
#' 
#' \dontrun{
#' list_cohorts()
#' }
#' @export
list_cohorts <- function() {
  cohort_info_url <- "http://gdac.broadinstitute.org/runs/info/firehose_get_disease_cohorts_list.txt"
  if (RCurl::url.exists(cohort_info_url)) {
    header = RCurl::dynCurlReader()
    stat_code <- RCurl::curlPerform(url = cohort_info_url, headerfunction = header$update, curl = header$curl())
    
    content <- header$value()[[1]]
    content <- strsplit(content, "\n")[[1]]
    
    idx <- grepl("^[A-Z]", content)
    ret <- content[idx]
    return(unlist(strsplit(ret, " ")))
  } else {
    stop(paste(cohort_info_url, "does not exist"), call. = FALSE)
  }
}

#' List available runs on firehose website.
#' 
#' List available runs on firehose website.
#' 
#' @return A character vector of runs.
#' @examples 
#' 
#' \dontrun{
#' list_runs()
#' }
#' @export
list_runs <- function() {
  run_info_url <- "http://gdac.broadinstitute.org/runs/info/firehose_get_public_runs_list.txt"
  if (RCurl::url.exists(run_info_url)) {
    header = RCurl::dynCurlReader()
    stat_code <- RCurl::curlPerform(url = run_info_url, headerfunction = header$update, curl = header$curl())
    
    content <- header$value()[[1]]
    content <- strsplit(content, "\n")[[1]]
    
    idx <- grepl("^[a-z]", content)
    ret <- content[idx]
    
    return(ret)
  } else {
    stop(paste(run_info_url, "does not exist"), call. = FALSE)
  }
}

#' List available run types on firehose website.
#' 
#' List available run types on firehose website.
#' 
#' @return A character vector of cohorts.
#' @examples 
#' 
#' \dontrun{
#' list_run_types()
#' }
#' @export
list_run_types <- function() {
  runs <- list_runs()
  ret <- unlist(strsplit(runs, "__"))
  
  ret <- ret[seq(1, length(ret), by = 2)]
  unique(ret)
}

#' Prepare mutation assessor dataset.
#' 
#' Prepare mutation assessor dataset.
#' 
#' @examples 
#' 
#' \dontrun{
#' prepare_ma()
#' }
#' @export
prepare_ma <- function() {
  ma_file <- paste0(path.package("bionexr"), "/madb.rds")
  if (file.exists(ma_file)) {
    return()
  }
  ma_url <- "https://s3-us-west-1.amazonaws.com/bionexr-data/madb.rds"
  download_file(ma_url)
  ret <- file.rename("madb.rds", ma_file)
}

download_file <- function(url, file_name = basename(url), resume = TRUE) {
  # if resume = TRUE and there is a local broken file, we resume download
  if (resume && file.exists(file_name)) {
    ret <- resume_download(url, broken_file = file_name)
    return(invisible(ret))
  }
  
  f = RCurl::CFILE(file_name, mode="wb")
  on.exit(RCurl::close(f))
  
  download_file_size <- get_download_file_size(url)
  pb <- txtProgressBar(1, download_file_size, style = 3)
  
  message("Downloading ", file_name, " -- total file size is ", human_readable(download_file_size), ", downloaded percentage:")
  ret = RCurl::curlPerform(url = url, 
                           writedata = f@ref, 
                           noprogress = FALSE, 
                           progressfunction = function(down, up) {
                             setTxtProgressBar(pb, down[[2]])
                           })
  invisible(list(stat_code = ret,
                 data_link = url,
                 file_name = file_name))
}

# Resume broken downloads.
# 
# Resume broken downloads.
# 
# @param url downloading url.
# @param broken_file broken file name.
# @examples 
# 
# \dontrun{
# resume_download(data_url)
# }
resume_download <- function(url, broken_file = basename(url)) {
  f = RCurl::CFILE(broken_file, mode="a+b")
  on.exit(RCurl::close(f))
  
  local_file_size <- get_local_file_size(broken_file)
  download_file_size <- get_download_file_size(url)
  
  if (local_file_size == download_file_size) {
    return(invisible(list(stat_code = 0,
                   data_link = url,
                   file_name = broken_file)))
  }
  
  pb <- txtProgressBar(1, download_file_size, style = 3)
  
  message("Resuming the download of ", broken_file, " -- total file size is ", human_readable(download_file_size), ", downloaded percentage:")
  ret = RCurl::curlPerform(url = url, 
                           writedata = f@ref, 
                           noprogress = FALSE, 
                           resume.from.large = local_file_size, 
                           progressfunction = function(down, up) {
                             setTxtProgressBar(pb, local_file_size + down[[2]])
                           })
  invisible(list(stat_code = ret,
                 data_link = url,
                 file_name = broken_file))
}

# Match data type
# 
# Match data type
match_data_type <- function(data_links, data_type) {
  file_names <- basename(data_links)
  
  # data_type in format: "run_type"_"data_type"
  switch(data_type, 
         # stddata =========================
         stddata_mutation = {
           pattern <- "Mutation_Packager_Calls\\.Level_3"
           idx <- grepl(pattern, file_names)
           if (sum(idx) != 1) {
             stop("Can not handle multiple match now", call. = FALSE)
           }
           
           data_links[idx]
         },
         
         stddata_expression = {
           pattern <- "Merge_rnaseqv2.+RSEM_genes__data\\.Level_3"
           idx <- grepl(pattern, file_names)
           if (sum(idx) == 0) {
             stop("Can not find expression file, please go to firehose to check", call. = FALSE)
           }
           if (sum(idx) > 2) {
             stop("Incorrect file number", call. = FALSE)
           }
           
           data_links[idx]
         },

         # analyses data ====================
         analyses_Aggregate_AnalysisFeatures = ,
         analyses = {
           pattern <- paste0(substring(data_type, 10), "\\.Level_4")
           idx <- grepl(pattern, file_names)
           if (sum(idx) != 1) {
             stop("Can not handle multiple match now", call. = FALSE)
           }
           data_links[idx]
         },
         # default action ====================
         {
           stop("Incorrect data type", call. = FALSE)
         })
}

# Get the file size on the url.
# 
# Get the file size on the url.
# 
# @param url A string of url.
# @examples 
# 
# \dontrun{
# get_download_file_size(data_url)
# }
get_download_file_size <- function(url) {
  h <- RCurl::basicHeaderGatherer()
  res <- RCurl::getURLContent(url, nobody = TRUE, headerfunction = h$update)
  ret <- h$value()[["Content-Length"]]
  
  as.integer(ret)
}
