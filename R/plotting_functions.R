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

#' Title
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
