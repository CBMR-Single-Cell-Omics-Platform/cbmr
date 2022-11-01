#' Remove non-breaking space
#'
#' @param x string
#'
#' @return string, where non-breaking space has been replaced with regular space
#' @export
fix_nonbreaking_space <- function(x) {
  gsub("(*UCP)\\s", " ", x, perl = TRUE)
}
