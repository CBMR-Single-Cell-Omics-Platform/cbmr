#' format numbers for exponent
#'
#' @return function that converts numbers to 10^x format
#' @export
#'
#' @examples
exponent_format <- function() {
  function(x) {
    parse(text=gsub("e\\+?", " %*% 10^", scales::scientific_format()(x)))
  }
}
