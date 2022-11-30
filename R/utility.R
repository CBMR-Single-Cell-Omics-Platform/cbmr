#' Remove non-breaking space
#'
#' @param x string
#'
#' @return string, where non-breaking space has been replaced with regular space
#' @export
fix_nonbreaking_space <- function(x) {
  gsub("(*UCP)\\s", " ", x, perl = TRUE)
}

#' Apply function line by line
#'
#' @param file path to file to be processed
#' @param fun function to be applied to line
#' @param ... other parameters applied to function
#'
#' @return
#' @export
#'
#' @examples
apply_by_line <- function(file, fun, ..., chunks = 10000) {
  connection <- file(file, open = "r")
  
  if (tools::file_ext(file) == "gz") {
    connection <- gzcon(connection)
  }
  
  out <- vector(mode = "list", length = chunks)
  n <- 1
  
  while(TRUE) {
    line <- readLines(con = connection, n = 1)
    if (length(line) == 0) {
      break
    }
    if (n %% chunks == 0) {
      out <- c(out, vector(mode = "list", length = chunks))
    }
    out[[n]] <- fun(line, ...)
    n <- n + 1
  }
  close(connection)
  out[lengths(out) != 0]
}
