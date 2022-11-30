#' Find nearest timepoint
#'
#' @param x vector, of new timepoints
#' @param timepoints vector, of target timepoints
#' @param cycle_length integer, length of a cycle. Typically 24 or 60.
#'
#' @return vector of same length as x, where the entries are the nearest value
#' supplied in 'timepoints'. In case of ties the first entry in 'timepoints'
#' is selected.
#' @export
#'
#' @examples
#' nearest_timepoint(0:36, timepoints = 0:5 * 4)
nearest_timepoint <- function(x, timepoints, cycle_length = 24) {
  rescale <- function(y) y * 2 * pi / cycle_length
  
  i <- complex(real = 0, imaginary = 1)
  
  timepoints_complex <- exp(rescale(timepoints) * i)
  x_complex <- exp(rescale(x) * i)
  
  distances <- Mod(outer(timepoints_complex, x_complex, FUN = `-`))
  # Floating point precision is finicky when two points are equally far apart.
  # Rounding leads to consistent handling.
  distances <- round(distances, digits = 8)
  
  idx_nearest <- apply(distances, MARGIN = 2, which.min)
  timepoints[idx_nearest]
}

#' Limorhyde function with multiple harmonics
#'
#' @param df data.frame 
#' @param time_colname 
#' @param harmonics 
#' @param period 
#' @param sinusoid 
#' @param n_knots 
#'
#' @return
#' @export
limorhyde_harmonics <- function(df, time_colname, harmonics = 1, period = 24, 
                                sinusoid = TRUE, n_knots = 3) 
{
  if (!(time_colname %in% colnames(df))) {
    stop("time_colname must be a named column in df.")
  }
  if (sinusoid) {
    get_vals <- function(harmonic, fun){
      fun(df[[time_colname]]/period * 2 * harmonic * pi)
    }
    cosvals <- lapply(seq_len(harmonics), get_vals, fun = cos)
    sinvals <- lapply(seq_len(harmonics), get_vals, fun = sin)
    d <- data.frame(cosvals, sinvals)
    new_colnames <- as.vector(outer(X = seq_len(harmonics), 
                         Y = c("_cos", "_sin"), 
                         FUN = "paste0"))
    colnames(d) <- paste(time_colname, new_colnames, sep = "_")
  }
  else {
    if (!require("bigsplines")) {
      stop("Package bigsplines required for fitting a spline. 
       Please install bigsplines and try again.")
    }
    knots <- seq(0, period - period/n_knots, length = n_knots)
    d <- bigsplines::ssBasis(df[[time_colname]] %% period, knots = knots, 
                            xmin = 0, xmax = period, periodic = TRUE)$X
    d <- as.data.frame(d)
    colnames(d) <- paste0(time_colname, "_knot", 1:n_knots)
  }
  return(d)
}