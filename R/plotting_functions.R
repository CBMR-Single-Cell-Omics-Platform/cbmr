#' format numbers for exponent
#'
#' @return function that converts numbers to 10^x format
#' @export
#'
#' @examples
#' exponent_format()
exponent_format <- function() {
  function(x) {
    parse(text=gsub("e\\+?", " %*% 10^", scales::scientific_format()(x)))
  }
}

#' Minus log10 transformation
#'
#' @param base 	a positive or complex number: the base with respect to which 
#' logarithms are computed. Defaults to 10.
#'
#' @return transformer, see scales::trans_new
#' @export
#'
#' @examples
#' minus_log_trans()
minus_log_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(
    name      = paste0("minus-log-", format(base)), 
    transform = trans, 
    inverse   = inv, 
    breaks    = scales::log_breaks(base = base), 
    domain    = c(1e-100, Inf)
    )
}

#' Find P-value cutoff
#'
#' @param x vector of p-values
#' @param cutoff desired FDR cutoff
#' 
#' @details 
#' This function takes a list of P-values, and calculates what the first P-value
#' above a certain FDR threshold would have had to be in order to have an FDR
#' less than said threshold. Uses an iterative approach to approximate, but can
#' be made arbitrarily precise.
#' Designed to be used when generating volcano-plots
#' 
#' @return
#' @export
#'
#' @examples
find_pvalue_cutoff <- function(x, cutoff = 0.05)
{
  if (any(x < 0)) stop("Negative P-values not allowed")
  pvals <- sort(x)
  fdr <- p.adjust(pvals, method = "fdr")
  
  if (!any(fdr < cutoff)){
    return(cutoff/(length(pvals)))
  }
  
  firstFailing <- max(which(fdr < cutoff)) + 1
  
  diffFn <- function(x)
  {
    p.adjust(c(x, pvals[-(firstFailing)]), method = "fdr")[[1]] - cutoff
  }
  optimize(diffFn, pvals[c(firstFailing - 1, firstFailing)])$minimum
}

#' Genomic log transformation
#' @details 
#' Returns a transformer that nicely handles log10 scale in a genomic context.
#' 
#' @return transformer, see scales::trans_new
#' @export
#'
#' @examples
genomic_log_trans <- function() {
  .mod_transform <- function(x) ifelse(x == 0, 0, sign(x) * (log10(abs(x)) + 1))
  .mod_inverse <- function(x) ifelse(x == 0, 0, sign(x) * 10^(abs(x) - 1))
  
  mod_breaks <- function(n = 8){
    function(x){
      breaks <- .mod_transform(x) %>%
        pretty(n = n) %>%
        .mod_inverse()
      return(breaks)
    }
  }
  
  labFn <- function(x)
  {
    out <- x
    x <- na.omit(x)
    suffix <- character(length(x))
    suffix[abs(x) < 10^3] <- "bp"
    suffix[abs(x) >= 10^3] <- "kb"
    suffix[abs(x) >= 10^6] <- "Mb"
    suffix[abs(x) >= 10^9] <- "Gb"
    
    ind1 <- abs(x) >= 10^3 & abs(x) < 10^6
    ind2 <- abs(x) >= 10^6 & abs(x) < 10^9
    ind3 <- abs(x) >= 10^9
    x[ind1] <- round(x[ind1]/10^3, 2)
    x[ind2] <- round(x[ind2]/10^6, 2)
    x[ind3] <- round(x[ind3]/10^9, 2)
    out[!is.na(out)] <- paste0(x, suffix)
    out
  }
  
  scales::trans_new(name = "genomic", 
                    inverse = .mod_inverse,
                    transform = .mod_transform,
                    breaks = mod_breaks(), 
                    format = labFn
  )
}

volcanoplot <- function(table = NULL, logfc_cutoff = NULL, fdr_cutoff = NULL,
                        format = "auto", logfc_col = NULL, pval_col = NULL, 
                        fdr_col = NULL) {
  if(is.null(table)) stop("table must be provided.")
  
  if (format == "auto") {
    if (all(c("logFC", "PValue", "FDR") %in% colnames(table))) {
      logfc_col <- "logFC"
      pval_col <- "PValue"
      fdr_col <- "FDR"
    } else if (all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(table))) {
      logfc_col <- "logFC"
      pval_col <- "P.Value"
      fdr_col <- "adj.P.Val"
    } else {
      stop("table type could not be automatically determined, ", 
      "use format = 'manual' instead.")
    }
  } else if (format == "manual") {
    if (is.null(logfc_col)) 
      stop("logfc_col must be specified when format is 'manual'")
    if(!(logfc_col %in% colnames(table))) {
      stop(logfc_col, " not found in table")
    }
    
    if (is.null(pval_col)) 
      stop("pval_col must be specified when format is 'manual'")
    if(!(pval_col %in% colnames(table))) {
      stop(pval_col, " not found in table")
    }
    
    if (!is.null(fdr_cutoff)){
      if (is.null(fdr_col))
        stop("fdr_col must be specified when format is 'manual' and ",
             "fdr_cutoff is specified")
      if(!(fdr_col %in% colnames(table))) {
        stop(fdr_col, " not found in table")
      }
    }
  } else {
    stop("format must be 'auto' or 'manual'.")
  }
  
  minLogFc <- min(table[[logfc_col]])
  maxLogFc <- max(table[[logfc_col]])
  minPval <- min(table[[pval_col]])
  maxPval <- 1
  
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = logfc_col, 
                      y = pval_col)) + 
    ggplot2::geom_point(size = 1, alpha = 0.5) +
    ggplot2::scale_y_continuous(name = "P-value", # expression(paste(-log[10], "(P-value)")
                                trans = minus_log_trans(), # Maybe we want to manually do the transformation
                                limits = c(maxPval, minPval)) +
    ggplot2::scale_x_continuous(name = expression(paste(log[2], 
                                                        "(Fold Change)")),
                                limits = c(minLogFc, maxLogFc)) +
    theme_bw()
  
  if (!is.null(logfc_cutoff) || !is.null(fdr_cutoff)) {
    if (!is.null(logfc_cutoff) && !is.null(fdr_cutoff)) {
      signif_string <- paste("abs(", logfc_col, ") >", logfc_cutoff, "&", 
                             fdr_col, "<", fdr_cutoff)
      
      legend_signif <- paste0("abs(", logfc_col, ") >", logfc_cutoff, " &\n",
                              "abs(", fdr_col, ") < ", fdr_cutoff)
      legend_nonsignif <- "Otherwise"
      p <- p + 
        ggplot2::geom_vline(xintercept = c(-1, 1) * logfc_cutoff) +
        ggplot2::geom_hline(
          yintercept = find_pvalue_cutoff(x      = table[[pval_col]], 
                                          cutoff = fdr_cutoff))
    } else if (!is.null(logfc_cutoff)) {
      signif_string <- paste0("abs(", logfc_col, ") > ", logfc_cutoff)
      legend_signif <- paste0("abs(", logfc_col, ") > ", logfc_cutoff)
      legend_nonsignif <- paste("abs(", logfc_col, ") ≤", logfc_cutoff)
      
      p <- p + 
        ggplot2::geom_vline(xintercept = c(-1, 1) * logfc_cutoff)
    } else {
      signif_string <- paste(fdr_col, "<", fdr_cutoff)
      legend_signif <- paste(fdr_col, "<", fdr_cutoff)
      legend_nonsignif <- paste(fdr_col, "≥", fdr_cutoff)
      
      p <- p + 
        ggplot2::geom_hline(
          yintercept = find_pvalue_cutoff(x      = table[[pval_col]], 
                                          cutoff = fdr_cutoff))
    }
    
    p <- p %+% ggplot2::aes_string(colour = signif_string)
    p <- p + ggplot2::scale_color_manual(values = c("TRUE" = "red", 
                                                    "FALSE" = "black"), 
                           name="Significance", 
                           labels = c("TRUE" = legend_signif, 
                                      "FALSE" = legend_nonsignif), 
                           breaks = c("TRUE", "FALSE"))
  }
  p
}
