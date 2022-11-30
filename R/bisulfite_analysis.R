#' Get methylation metric
#'
#' @description 
#' Get various per sample methylation metrics. Currently only percent methylation
#' and Mvalue (methylated counts / unmethylated counts) are supported.
#'
#' @param x object to calculate methylation metric from
#' @param metric string, method to use, currently only percent and Mvalue are supported
#' @export
#'
#' @details 
#' x can be a matrix, a data.frame, a data.table or an edgeR DGEList.
#' 
#' x must contain only columns with counts and each sample must have two columns,
#' one ending in '_Me' and one ending in '_Un' with methylated and unmethylated
#' counts respectively. The columns need not be in that order.
#'
#' @return object of the same type as x, except for DGEList which returns an object of
#' the same type as x$counts (typically a matrix).
#' @examples 
#' set.seed(1)
#' n_samples <- 4
#' n_genes <- 10
#' counts <- matrix(sample(1:20, n_samples * 2 * n_genes, replace = TRUE), 
#'                  ncol = n_samples * 2)
#' colnames(counts) <- paste0("sample", rep(1:4, each = 2), c("_Me", "_Un"))
#' get_methylation_metric(counts)
#' get_methylation_metric(counts, metric = "Mvalue")
get_methylation_metric <- function(x, metric = "percent") {
  UseMethod("get_methylation_metric")
}

#' @describeIn get_methylation_metric Methylation metric default
#' @export
get_methylation_metric.default <- function(x, metric = "percent"){
  sample_names <- unique(stringr::str_remove(colnames(x), "_Me$|_Un$"))
  
  missed_cols <- setdiff(paste0(rep(sample_names, each = 2), c("_Me", "_Un")), colnames(x))
  if (length(missed_cols) != 0) stop('Counts matrix must contain two columns for each sample, one ending in "_Me" and one ending in "_Un"')
  
  me_idx <- paste0(sample_names, "_Me")
  un_idx <- paste0(sample_names, "_Un")
  
  if (metric == "percent") {
    out <- x[, me_idx] / (x[, me_idx] + x[, un_idx])
  } else if (metric == "Mvalue") {
    out <- x[, me_idx] / x[, un_idx]
  } else {
    stop('metric must be one of "percent" and "Mvalue"')
  }
  colnames(out) <- sample_names
  out
}

#' @describeIn get_methylation_metric Methylation metric for data.tables
#' @export
get_methylation_metric.data.table <- function(x, metric = "percent"){
  sample_names <- unique(stringr::str_remove(colnames(x), "_Me$|_Un$"))
  
  missed_cols <- setdiff(paste0(rep(sample_names, each = 2), c("_Me", "_Un")), colnames(x))
  if (length(missed_cols) != 0) stop('Counts matrix must contain two columns for each sample, one ending in "_Me" and one ending in "_Un"')
  
  me_idx <- paste0(sample_names, "_Me")
  un_idx <- paste0(sample_names, "_Un")
  
  if (metric == "percent") {
    out <- x[, me_idx, with = FALSE] / (x[, me_idx, with = FALSE] + x[, un_idx, with = FALSE])
  } else if (metric == "Mvalue") {
    out <- x[, me_idx, with = FALSE] / x[, un_idx, with = FALSE]
  } else {
    stop('metric must be one of "percent" and "Mvalue"')
  }
  data.table::setnames(out, sample_names)
  out
}

#' @describeIn get_methylation_metric Methylation metric for DGELists
#' @export
get_methylation_metric.DGEList <- function(x, metric = "percent"){
  get_methylation_metric(x$counts, metric = metric)
}

#' Get groupwise average methylation metric
#'
#' @description Get the average groupwise methylation metric (in case of a '~ 0 + group' desing)
#' or effect size (in case of a '~ a * b' model).
#'
#' @param methylation matrix or object coercible to a matrix containing methylation metrics
#' @param design design matrix with experimental design
#'
#' @return matrix with groupwise averages of methylation metrics
#' @export
#'
#' @examples
#' set.seed(1)
#' n_samples <- 4
#' n_genes <- 10
#' counts <- matrix(sample(1:20, n_samples * 2 * n_genes, replace = TRUE), 
#'                  ncol = n_samples * 2)
#' colnames(counts) <- paste0("sample", rep(1:4, each = 2), c("_Me", "_Un"))
#' methylation_levels <- get_methylation_metric(counts)
#' 
#' group <- c("a", "a", "b", "b")
#' design <- model.matrix(~ 0 + group)
#' 
#' get_groupwise_methylation(methylation_levels, design)
get_groupwise_methylation <- function(methylation, design) {
  if (!inherits(methylation, "matrix")){
    warning("Coercing methylation to matrix")
    methylation <- as.matrix(methylation)
  }
  # Standard estimator of beta_hat from linear regression
  methylation %*% t(solve(t(design) %*% design) %*% t(design))
}

#' Get delta methylation metric
#'
#' @description Calculate effect sizes of contrasts
#'
#' @param group_methylation object coercible to a matrix containing groupwise 
#' methylation (or effect sizes) as calcualted by get_groupwise_methylation
#' @param contrasts matrix of contrasts as created by edgeR::makeContrasts
#'
#' @return matrix with the effect sizes of contrasts
#' @export
#'
#' @examples
#' require(limma)
#' set.seed(1)
#' n_samples <- 4
#' n_genes <- 10
#' counts <- matrix(sample(1:20, n_samples * 2 * n_genes, replace = TRUE), 
#'                  ncol = n_samples * 2)
#' colnames(counts) <- paste0("sample", rep(1:4, each = 2), c("_Me", "_Un"))
#' methylation_levels <- get_methylation_metric(counts)
#' 
#' group <- c("a", "a", "b", "b")
#' design <- model.matrix(~ 0 + group)
#' 
#' group_methylation <- get_groupwise_methylation(methylation_levels, design)
#' 
#' contrasts <- makeContrasts(groupb - groupa, (groupb + groupa)/2, levels = design)
#' get_delta_methylation(group_methylation, contrasts)

get_delta_methylation <- function(group_methylation, contrasts) {
  if (!inherits(group_methylation, "matrix")){
    warning("Coercing group_methylation to matrix")
    group_methylation <- as.matrix(group_methylation)
  }
  
  same_names <- all(colnames(group_methylation) == rownames(contrasts))
  if (!same_names) {
    colnames(group_methylation)[colnames(group_methylation) == "(Intercept)"] <- "Intercept"
    same_names <- all(colnames(group_methylation) == rownames(contrasts))
    if (!same_names) {
      stop("Column names of group methylation must match rownames of contrast")
    } else {
      warning("Renaming (Intercept) to Intercept")
    }
  }
  
  group_methylation %*% contrasts
}


#' Count filtering for methylation counts
#'
#' @param y DGEList
#' @param min_counts minimum number of counts
#' @param filter_extremes logical, filter 0% or 100% in all
#'
#' @return filtered DGEList
#' @export
meth_filter <- function(y, min_counts = 8, filter_extremes = TRUE) {
  UseMethod("meth_filter")
}

#' @describeIn meth_filter Count filtering for methylation counts DGEList
#' @export
meth_filter.DGEList <- function(y, min_counts = 8, filter_extremes = TRUE)
{
  me_idx <- stringr::str_subset(colnames(y), "_Me$")
  un_idx <- stringr::str_subset(colnames(y), "_Un$")
  
  # Filter by counts
  coverage <- y$counts[, me_idx] +
    y$counts[, un_idx]
  
  keep <- rowSums(coverage >= min_counts) == ncol(coverage)
  y <- y[keep,, keep.lib.sizes=FALSE]
  
  if (filter_extremes) {
    # From discussion with Gordon, this should be OK
    constant <- rowSums(y$counts[, me_idx]) == 0 |
      rowSums(y$counts[, un_idx]) == 0
    
    y <- y[!constant,, keep.lib.sizes=FALSE]
  }
  y
}

#' @describeIn meth_filter Count filtering for methylation counts data.frame
#' @export
meth_filter.data.frame <- function(y, min_counts = 8, filter_extremes = TRUE)
{
  me_idx <- stringr::str_detect(colnames(y), "_Me$")
  un_idx <- stringr::str_detect(colnames(y), "_Un$")
  
  # Filter by counts
  coverage <- y[, me_idx] +
    y[, un_idx]
  
  keep <- rowSums(coverage >= min_counts) == ncol(coverage)
  y <- y[keep,]
  
  if (filter_extremes) {
    # From discussion with Gordon, this should be OK
    constant <- rowSums(y[, (me_idx)]) == 0 |
      rowSums(y[, (un_idx)]) == 0
    
    y <- y[!constant,,]
  }
  y
}

#' @describeIn meth_filter Count filtering for methylation counts data.table
#' @export
meth_filter.data.table <- function(y, min_counts = 8, filter_extremes = TRUE)
{
  # Filter by counts
  coverage <- y[, .SD, .SDcols = patterns("_Me$")] +
    y[, .SD, .SDcols = patterns("_Un$")]
  
  keep <- rowSums(coverage >= min_counts) == ncol(coverage)
  y <- y[keep,]
  
  if (filter_extremes) {
    # From discussion with Gordon, this should be OK
    constant <- rowSums(y[, .SD, .SDcols = patterns("_Me$")]) == 0 |
      rowSums(y[, .SD, .SDcols = patterns("_Un$")]) == 0
    
    y <- y[!constant,,]
  }
  y
}

#' Wrapper for running RRBS analysis through edgeR
#'
#' @param y DGEList with a column for unmethylated and a column for methylated counts per sample.
#' @param design design as generated by model.matrix. modelMatrixMeth is called in this function.
#' @param parallel logical, should large datasets be run in parallel. If true, samples with more than 
#' 100.000 rows are split into 50.000 row chunks and analyzed separately.
#' @param aggregate_by optional, name of column in y$genes to aggregate by.
#' @param sample_pattern 
#'
#' @return
#' @export
#' @importFrom foreach %dopar%
edger_for_rrbs <- function(y, design, parallel = TRUE)
{
  design <- modelMatrixMeth(design)
  
  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(floor(cores[1]/2))
  doParallel::registerDoParallel(cl)
  
  if(nrow(y) > 99999){
    len_out <- floor(nrow(y) / 50000)
    idxs <- ceiling(seq(1, nrow(y), length.out = len_out))
    idxs <- cbind(idxs[-length(idxs)] + 1, idxs[-1])
    idxs[1,1] <- idxs[1,1] - 1
  } else {
    idxs <- matrix(c(1, nrow(y)), nrow = 1)
  }
  
  y <- foreach::foreach(n = seq_len(nrow(idxs)), .packages = "edgeR") %dopar% {
    edgeR::estimateDisp(y[seq(idxs[n, 1], idxs[n, 2]), ], design=design, trend="none")
  }
  
  fit <- foreach::foreach(n = y, .packages = "edgeR") %dopar% {
    edgeR::glmFit(n, design)
  }
  
  parallel::stopCluster(cl)
  
  fit
}