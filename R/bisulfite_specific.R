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
  
  ..me_idx <- ..un_idx <- NULL
  
  if (metric == "percent") {
    out <- x[, ..me_idx] / (x[, ..me_idx] + x[, ..un_idx])
  } else if (metric == "Mvalue") {
    out <- x[, ..me_idx] / x[, ..un_idx]
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
  # Stanard estimator of beta_hat from linear regression
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
#' contrasts <- limma::makeContrasts(groupb - groupa, (groupb + groupa)/2, levels = design)
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
