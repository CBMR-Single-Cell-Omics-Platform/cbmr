#' Prepare density plots
#'
#' @param y DGEList
#'
#' @return One or more plots are created on the current device, 
#' list of the plot generated are also returned
#' @export
#'
#' @examples
#' library(edgeR)
#' set.seed(1)
#' counts <- matrix(sample(1:20), ncol = 4)
#' y <- DGEList(counts)
#' p <- plot_density(y)
#' p
plot_density <- function(y)
{
  pD <- reshape2::melt(edgeR::cpm(y, normalized.lib.sizes = TRUE, log = TRUE), 
                       varname = c("Var1", "Sample"), as.is = TRUE)
  pD$Sample <- factor(pD$Sample, levels = colnames(y))
  
  if (ncol(y) > 16) {
    
    sample_groups <- split(colnames(y), ceiling(seq_along(colnames(y))/16))
    sample_groups <- utils::stack(sample_groups)
    sample_groups <- structure(sample_groups$ind, names = sample_groups$values)
    
    samples <- split(pD, sample_groups[pD$Sample])
    
    value <- NULL
    
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
#' @param y DGEList
#'
#' @return One or more plots are created on the current device
#' @export
#'
#' @examples
#' library(edgeR)
#' set.seed(1)
#' counts <- matrix(sample(1:20), ncol = 4)
#' y <- DGEList(counts)
#' plot_all_md(y)
plot_all_md <- function(y)
{
  withr::with_par(list(mfrow = c(2, 2)), {
    for (i in seq_len(ncol(y))){
      edgeR::plotMD.DGEList(y, column = i)
      graphics::abline(h = 0)
    }
  })
}


#' Prepare data for MDS plot
#'
#' @param y DGEList, EList or matrix.
#' @param dim_plot integer vector of length two specifying which principal components should be plotted.
#' @param colour_by a column in y$samples or a vector of the same length as ncol(y). Set to NULL to disable colours.
#' @param metadata optional data.frame with additional meta data.
#'
#' @return data.frame with data necessary for an MDS plot
prepare_mds_data <- function(y, dim_plot, colour_by = NULL, metadata = NULL) {
  UseMethod("prepare_mds_data")
}

#' @describeIn prepare_mds_data Prepare data for MDS plot DGEList
#' @export
prepare_mds_data.EList <- function(y, dim_plot, colour_by = NULL, metadata = NULL) {
  plot_data <- limma::plotMDS(y, dim.plot = dim_plot, plot = FALSE)$eigen.vectors[, dim_plot]
  plot_data <- as.data.frame(plot_data)
  colnames(plot_data) <- paste0("dim", dim_plot)
  plot_data <- cbind(plot_data, y$targets)
  
  if (!is.null(metadata)) {
    plot_data <- cbind(plot_data, metadata)
  }
  
  # Colour by is either a column in y$samples or a vector of group memberships
  if (is.null(colour_by)) {
    colour_idx <- NULL
  } else if (length(colour_by) == 1) {
    if (colour_by %in% colnames(plot_data)) {
      colnames(plot_data)[colnames(plot_data) == colour_by] <- "colour_column"
    } else {
      stop(paste(colour_by, "not found in y$targets"))
    }
  } else if (length(colour_by) == ncol(y)) {
    plot_data <- cbind(plot_data, colour_column = colour_by)
  } else {
    stop("colour_by must be NULL, a column in y$targets or a vector of the same length as ncol(y)")
  }
  plot_data
}


#' @describeIn prepare_mds_data Prepare data for MDS plot DGEList
#' @export
prepare_mds_data.DGEList <- function(y, dim_plot, colour_by = NULL, metadata = NULL) {
  plot_data <- edgeR::plotMDS.DGEList(y, dim.plot = dim_plot, plot = FALSE)$eigen.vectors[, dim_plot]
  plot_data <- as.data.frame(plot_data)
  colnames(plot_data) <- paste0("dim", dim_plot)
  plot_data <- cbind(plot_data, y$samples)
  
  if (!is.null(metadata)) {
    plot_data <- cbind(plot_data, metadata)
  }
  
  # Colour by is either a column in y$samples or a vector of group memberships
  if (is.null(colour_by)) {
    colour_idx <- NULL
  } else if (length(colour_by) == 1) {
    if (colour_by %in% colnames(plot_data)) {
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

#' @describeIn prepare_mds_data Prepare data for MDS plot default
#' @export
prepare_mds_data.default <- function(y, dim_plot, colour_by = NULL, metadata = NULL) {
  plot_data <- limma::plotMDS(y, dim.plot = dim_plot, plot = FALSE)$eigen.vectors[, dim_plot]
  plot_data <- as.data.frame(plot_data)
  colnames(plot_data) <- paste0("dim", dim_plot)
  
  if (!is.null(metadata)) {
    plot_data <- cbind(plot_data, metadata)
  }
  
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
#' @param y DGEList, EList or object that can be handled by limma::plotMDS
#' @param dim_plot integer vector of length two specifying which principal components should be plotted.
#' @param colour_by a column in y$samples or a vector of the same length as ncol(y). Set to NULL to disable colours.
#' @param col_scale optional colour scale to be used when plotting.
#' @param size numeric, size of points
#'
#' @importFrom ggplot2 %+%
#' 
#' @description 
#' Prepares a ggplot2 object with the data obtained from running limma/edgeR
#' plotMDS. If y is a DGEList the data in y$samples is appended so the plot can be 
#' updated to reflect different columns in the metadata. See examples.
#' 
#' @return ggplot2 object which when evaluated creates a MDS plot.
#' @export
#'
#' @examples
#' library(edgeR)
#' library(ggplot2)
#' set.seed(1)
#' counts <- matrix(sample(1:20), ncol = 4)
#' samples <- cbind(
#'   foo = letters[c(1,1,2,2)],
#'   bar = LETTERS[c(1,2,1,2)]
#' )
#' y <- DGEList(counts, samples = samples)
#' p <- ggplot_mds(y, dim_plot = 1:2)
#' 
#' p
#' p %+% aes(colour = foo)
#' p %+% aes(colour = bar)

ggplot_mds <- function(y, dim_plot, colour_by = NULL, col_scale, size = 5) {
  
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
    ggplot2::geom_point(size = size) +
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

#' Prepare data for heatmap
#'
#' @param x object containing read counts, either a matrix, a DGEList object or an EList object
#' @param method character, method to calculate distances, one of 'MDS', 
#' 'poission' and the methods supported by stats::dist, see details.
#' @param cpm logical, converted to counts per million? 
#' Only used for the methods supported by stats::dist
#' @param log logical, log transform?
#' Only used for the methods supported by stats::dist
#'
#' @details 
#' Uses a number of methods to prepare data for creating a heatmap. x should 
#' contain read counts. If the method 'MDS' or 'poisson' is selected the raw counts
#' are processed with limma::plotMDS or PoiClaClu::PoissonDistance, otherwise counts
#' are optionally converted to counts per million and/or log transformed before
#' distances are calculated with stats::dist.
#'
#' @return a 'dist' object with sample to sample distances
#' @export
#'
#' @examples
#' library(edgeR)
#' set.seed(1)
#' counts <- matrix(sample(1:20), ncol = 4)
#' prepare_heatmap_data(counts, method = "euclidean")
prepare_heatmap_data <- function(x, method, cpm = FALSE, log = FALSE) {
  UseMethod("prepare_heatmap_data")
}

#' @describeIn prepare_heatmap_data prepare heatmap data default
#' @export
prepare_heatmap_data.default <- function(x, method, cpm = FALSE, log = FALSE) {
  if (method == "MDS") {
    out <- stats::as.dist(limma::plotMDS.default(x, dim.plot = c(1, 2), 
                                                     plot = FALSE)$distance.matrix)
  } else if (method == "poisson") {
    out <- PoiClaClu::PoissonDistance(t(x))$dd
    attr(out, "Labels") <- colnames(x)
  } else {
    if (cpm) {
      x <- t(edgeR::cpm(x, log = log))
    } else {
      x <- t(x)
      if (log) {
        if (any(x <= 0)) {
          stop("data contains 0 or negative values, log transform not supported")
        }
        x <- log(x)
      }
    }
    out <- stats::dist(x, method = method)
  }
  out
}

#' @describeIn prepare_heatmap_data prepare heatmap data for DGEList
#' @export
prepare_heatmap_data.DGEList <- function(x, method, cpm = FALSE, log = FALSE) {
  # It is tempting to just extract x$counts and pass to prepare_heatmap_data.default
  # but plotMDS and cpm can use the extra information encoded in the DGEList so 
  # we have to have a bit of code duplication
  if (method == "MDS") {
    out <- stats::as.dist(edgeR::plotMDS.DGEList(x, dim.plot = c(1, 2), 
                                                     plot = FALSE)$distance.matrix)
  } else if (method == "poisson") {
    out <- PoiClaClu::PoissonDistance(t(x$counts))$dd
    attr(out, "Labels") <- colnames(x)
  } else {
    if (cpm) {
      x <- t(edgeR::cpm(x, log = log))
    } else {
      x <- t(x$counts)
      if (log) {
        if (any(x <= 0)) {
          stop("data contains 0 or negative values, log transform not supported")
        }
        x <- log(x)
      }
    }
    out <- stats::dist(x, method = method)
  }
  out
}

#' @describeIn prepare_heatmap_data prepare heatmap data for DGEList
#' @export
prepare_heatmap_data.EList <- function(x, method, cpm = FALSE, log = FALSE) {
  prepare_heatmap_data.default(x$E, method, cpm = FALSE, log = FALSE)
}


#' Plot heatmap of sample distances
#'
#' @param x a DGE list or a matrix of numbers
#' @param method string, method used for calculating sample distances. Supported methods
#' are poisson, MDS, and the methods supported by dist. See details.
#' @param cpm logical, should input be converted to counts per million, see details
#' @param log logical, should CPM be log-transformed, see details
#' @param ... passed to pheatmap
#'
#' @details If the method is MDS, the limma/edgeR plotMDS is used to calculate
#' sample distances using the first two dimensions. If the method is poisson distance is calculated using the
#' poiClaClu package. Otherwise input is optionally converted to (log)CPM and 
#' passed to dist with the specified method.
#' 
#' If MDS or poisson is used, read counts should be supplied. 
#' 
#' @return pheatmap object
#' @export
#'
#' @examples
#' library(edgeR)
#' set.seed(1)
#' counts <- matrix(sample(1:20), ncol = 4)
#' plot_data <- plot_sample_heatmap(counts, method = "euclidean")
plot_sample_heatmap <- function(x, method, cpm = FALSE, log = FALSE, ...) {
  heatmap_data <- prepare_heatmap_data(x, method = method, cpm, log)
  
  pheatmap::pheatmap(as.matrix(heatmap_data), 
           clustering_distance_rows=heatmap_data,
           clustering_distance_cols=heatmap_data,
           legend = FALSE,
           ...)
}

#' Plot number of CpGs retained by including samples.
#'
#' @param x list of data.tables with CpG information
#'
#' @return ggplot object that creates the plot
#' @export
#' @examples 
#' plot_retained_cpgs(cov_list)

plot_retained_cpgs <- function(x) {
  n_cpg <- unlist(lapply(x, nrow))
  if (is.null(names(x))) names(x) <- paste("Sample", seq_along(x))
  idx <- order(n_cpg, decreasing = TRUE) #Assume cutoff is made based on number of CpGs covered
  
  n_cpg <- n_cpg[idx]
  x <- x[idx]
  
  
  merger <- function(x, y) data.table::merge.data.table(x, y, by = c("seqnames", "start", "end"))
  . <- seqnames <- start <- end <- sample_name <- NULL
  stepwise_retained <- Reduce(merger, lapply(x, `[`, , .(seqnames, start, end)), accumulate = TRUE)
  n_retained <- unlist(lapply(stepwise_retained, nrow))
  
  pD <- data.frame(n_cpg = n_cpg,
                   n_retained = n_retained, 
                   sample_name = names(x)
  )
  
  suppressWarnings(ggplot2::ggplot(pD, ggplot2::aes(x = n_cpg/10^6, y = n_retained/10^6, label = sample_name)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel() +
    ggplot2::scale_x_log10(name = "CpGs covered in sample [millions]") +
    ggplot2::scale_y_log10(name = "CpGs covered in all samples [millions]") +
    ggplot2::theme_bw()
  )
}
