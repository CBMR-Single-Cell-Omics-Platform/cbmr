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
    format    = scales::label_math(expr = 10^-.x, format = trans),
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

#' Volcanoplot
#'
#' @param table data.frame with output from edgeR or limma. Must contain three 
#' columns named "logFC", "PValue"/"P.Value" & "FDR"/"adj.P.Val"
#' @param logfc_cutoff numeric, absolute logFC must be higher than this to be
#' coloured significant
#' @param fdr_cutoff numeric, FDR/adj.P.Val must be less than this to be 
#' coloured significant. Uses \link[cbmr]{find_pvalue_cutoff} to determine where
#' to place cutoff line.
#' @param extra_pval numeric, P-value to include when calculating y-limits. 
#' Useful for setting the same limits in multiple volcanoplots
#' @param extra_logfc numeric, logFC to include when calculating x-limits. 
#' Useful for setting the same limits in multiple volcanoplots
#'
#' @importFrom ggplot2 %+%
#'
#' @export
volcanoplot <- function(table = NULL, logfc_cutoff = NULL, fdr_cutoff = NULL,
                        extra_pval = NULL, extra_logfc = NULL) {
  if(is.null(table)) stop("table must be provided.")
  
  if (all(c("logFC", "PValue", "FDR") %in% colnames(table))) {
    logfc_col <- "logFC"
    pval_col <- "PValue"
    fdr_col <- "FDR"
  } else if (all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(table))) {
    logfc_col <- "logFC"
    pval_col <- "P.Value"
    fdr_col <- "adj.P.Val"
  } else {
    stop("table type could not be automatically determined.")
  }
  
  minLogFc <- min(c(table[[logfc_col]], extra_logfc))
  maxLogFc <- max(c(table[[logfc_col]], extra_logfc))
  minPval <- min(c(table[[pval_col]], extra_pval))
  
  maxPval <- 1
  
  p <- ggplot2::ggplot(table, ggplot2::aes_string(x = logfc_col, 
                      y = pval_col)) + 
    ggplot2::geom_point(size = 1, alpha = 0.5) +
    ggplot2::scale_y_continuous(name = "P-value", # expression(paste(-log[10], "(P-value)")
                                trans = minus_log_trans(), # Maybe we want to manually do the transformation
                                limits = c(maxPval, minPval)) +
    ggplot2::scale_x_continuous(name = expression(paste(log[2], 
                                                        "(Fold Change)")),
                                limits = c(minLogFc, maxLogFc)) +
    ggplot2::theme_bw()
  
  if (!is.null(logfc_cutoff) || !is.null(fdr_cutoff)) {
    signif_string <- NULL
    nonsig_string <- NULL
    
    if (!is.null(logfc_cutoff)) {
      signif_string <- paste0("abs(", logfc_col, ") > ", logfc_cutoff)
      nonsig_string <- paste0("abs(", logfc_col, ") ≤ ", logfc_cutoff)
      
      p <- p + 
        ggplot2::geom_vline(xintercept = c(-1, 1) * logfc_cutoff)
    }
    
    if (!is.null(fdr_cutoff)) {
      signif_string <- c(signif_string, paste0(fdr_col, " < ", fdr_cutoff))
      nonsig_string <- c(nonsig_string, paste0(fdr_col, " ≥ ", fdr_cutoff))
      
      p <- p + 
        ggplot2::geom_hline(
          yintercept = find_pvalue_cutoff(x      = table[[pval_col]], 
                                          cutoff = fdr_cutoff))
    }
    
    legend_signif <- paste(signif_string, collapse = " &\n")
    legend_nonsig <- ifelse(test = length(nonsig_string) > 1, 
                            yes = "Otherwise", 
                            no = nonsig_string)
    
    signif_string <- paste(signif_string, collapse = " & ")
    
    p <- p %+% ggplot2::aes_string(colour = signif_string)
    p <- p + ggplot2::scale_color_manual(
      values = c("TRUE" = "red", 
                 "FALSE" = "black"), 
      name   = "Significance", 
      labels = c("TRUE"  = legend_signif, 
                 "FALSE" = legend_nonsig), 
      breaks = c("TRUE", "FALSE"))
  }
  p
}
