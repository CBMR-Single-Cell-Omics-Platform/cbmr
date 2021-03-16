#' Wrapper for testing coefficients or contrasts from an edgeR model
#'
#' @param contrast character, name of contrast to be tested
#' @param coef character, name of coefficient to be tested
#' @param efit DGEGLM object, object as created by edgeR function glmQLFit or glmFit
#' @param contast_matrix matrix, typically created by the edgeR function makeContrasts
#' @param id_name character, name of the gene ids
#'
#' @return
#' @import data.table
#' @export
#'
#' @examples
edgeR_tester <- function(contrast = NULL, coef = NULL, efit, contast_matrix, id_name = "ENSEMBL") {
  if (missing(efit)) stop("efit must be specified")
  
  if (!is.null(coef)) {
    qlfTest <- edgeR::glmQLFTest(glmfit = efit, coef = coef)
  } else if (!is.null(contrast)) {
    if (missing(contast_matrix)) stop("contast_matrix must be specified when contrast is used")
    qlfTest <- edgeR::glmQLFTest(glmfit = efit, contrast = contast_matrix[, contrast])
  } else {
    stop("Must specify either coef or contrast")
  }
  out <- edgeR::topTags(qlfTest, n = Inf, p.value = 1)$table
  data.table::setDT(out, keep.rownames = TRUE)
  data.table::setnames(out, "rn", id_name)
  out[]
}
