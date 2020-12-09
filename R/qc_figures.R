#' Prepare density plots
#'
#' @param y DGEList to be plotted
#'
#' @return One or more plots are created on the current device, 
#' list of the plot generated are also returned
#' @export
#'
#' @examples
#' \dontrun{
#'   p <- plot_density(y)
#'   p
#' }
plot_density <- function(y)
{
  pD <- reshape2::melt(edgeR::cpm(y, normalized.lib.sizes = TRUE, log = TRUE), 
                       varname = c("Var1", "Sample"), as.is = TRUE)
  pD$Sample <- factor(pD$Sample, levels = colnames(y))
  
  if (ncol(y) > 16) {
    
    sample_groups <- split(colnames(y), ceiling(seq_along(colnames(y))/16))
    sample_groups <- stack(sample_groups)
    sample_groups <- structure(sample_groups$ind, names = sample_groups$values)
    
    samples <- split(pD, sample_groups[pD$Sample])
    
    p <- vector(mode = "list", length = length(samples))
    for (i in seq_along(samples)) {
      p[[i]] <- ggplot2::ggplot(samples[[i]], ggplot2::aes(value)) + 
        ggplot2::geom_density() + 
        ggplot2::facet_wrap(~Sample)
    }
  } else {
    p <- ggplot2::ggplot(pD, ggplot2::aes(value)) + 
      ggplot2::geom_density() + 
      ggplot2::facet_wrap(~Sample)
  }
  invisible(p)
}

#' Prepare MD plots
#'
#' @param y 
#'
#' @return One or more plots are created on the current device
#' @export
#'
#' @examples
#' \dontrun{
#'   plot_all_md(y)
#' }
plot_all_md <- function(y)
{
  withr::with_par(list(mfrow = c(2, 2)), {
    for (i in seq_len(ncol(y))){
      edgeR::plotMD.DGEList(y, column = i)
      graphics::abline(h = 0)
    }
  })
}

prepare_mds_data <- function(y,...) UseMethod("prepare_mds_data")

#' Prepare data for MDS plot
#'
#' @param y 
#' @param dim_plot 
#' @param colour_by 
#' @param col_scale 
#'
#' @return
#' @export
#'
#' @examples
prepare_mds_data.DGEList <- function(y, dim_plot, colour_by = NULL) {
  plot_data <- edgeR::plotMDS.DGEList(y, dim.plot = dim_plot, plot = FALSE)$cmdscale.out
  plot_data <- as.data.frame(plot_data)
  colnames(plot_data) <- paste0("dim", dim_plot)
  
  # Colour by is either a column in y$samples or a vector of group memberships
  if (is.null(colour_by)) {
    colour_idx <- NULL
  } else if (length(colour_by) == 1) {
    if (colour_by %in% colnames(y$samples)) {
      plot_data <- cbind(plot_data, y$samples)
      colnames(plot_data)[colnames(plot_data) == colour_by] <- "colour_column"
    } else {
      stop(paste(colour_by, "not found in y$samples"))
    }
  } else if (length(colour_by) == ncol(y)) {
    plot_data <- cbind(plot_data, colour_column = colour_by)
  } else {
    stop("colour_by must be NULL, a column in y$samples or a vector of the same length as ncol(y)")
  }
  plot_data
}

#' Prepare data for MDS plot
#'
#' @param y 
#' @param dim_plot 
#' @param colour_by 
#' @param col_scale 
#'
#' @return
#'
#' @examples
prepare_mds_data.default <- function(y, dim_plot, colour_by = NULL) {
  plot_data <- limma::plotMDS.default(y, dim.plot = dim_plot, plot = FALSE)$cmdscale.out
  plot_data <- as.data.frame(plot_data)
  colnames(plot_data) <- paste0("dim", dim_plot)
  
  # Colour by is either a column in y$samples or a vector of group memberships
  if (is.null(colour_by)) {
    colour_idx <- NULL
  } else if (length(colour_by) == ncol(y)) {
    plot_data <- cbind(plot_data, colour_column = colour_by)
  } else {
    stop("colour_by must be NULL or a vector of the same length as ncol(y)")
  }
  plot_data
}

#' Plot MDS plot
#'
#' @param y 
#' @param dim_plot 
#' @param colour_by 
#' @param col_scale 
#'
#' @importFrom ggplot2 %+%
#'
#' @return
#' @export
#'
#' @examples
ggplot_mds <- function(y, dim_plot, colour_by = NULL, col_scale) {
  
  plot_data <- prepare_mds_data(y, dim_plot, colour_by = colour_by)
  if (!is.null(colour_by)) {
    colour <- "colour_column"
  } else {
    colour <- NULL
  }
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = colnames(plot_data)[1],
                                                      y = colnames(plot_data)[2],
                                                      colour = colour)
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_bw()
  
  if (!is.null(colour_by)){
    if (!missing(col_scale)) {
      p <- p + col_scale
    } else if (length(colour_by) == 1){
      p <- p %+% ggplot2::labs(colour = colour_by)
    } else {
      p <- p %+% ggplot2::labs(colour = NULL)
    }
  }
  p
}