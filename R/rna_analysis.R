#' Wrapper for testing coefficients or contrasts from an edgeR model
#'
#' @param contrast character, name of contrast to be tested
#' @param coef character, name of coefficient to be tested
#' @param efit DGEGLM object, object as created by edgeR function glmQLFit or glmFit
#' @param contrast_matrix matrix, typically created by the edgeR function makeContrasts
#' @param id_name character, name of the gene ids
#'
#' @return data.table with edgeR results
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' y <- prepare_rna_dgelist(rna_counts, example_metadata, sample_col = "Sample.ID")
#' design <- model.matrix(~ 0 + Genotype + Batch, data = y$samples)
#' y <- edgeR::estimateDisp(y, design = design)
#' rownames(y) <- y$genes$Geneid
#' contrast_matrix <- limma::makeContrasts(Genotype = GenotypeKO - GenotypeWT, 
#'                                         levels = y$design)
#' efit <- glmQLFit(y, design = y$design, robust = TRUE)
#' edgeR_tester("Genotype", efit = efit, contrast_matrix = contrast_matrix)
#' }
edgeR_tester <- function(contrast = NULL, coef = NULL, efit, contrast_matrix, id_name = "ENSEMBL") {
  if (missing(efit)) stop("efit must be specified")
  
  if (!is.null(coef)) {
    qlfTest <- edgeR::glmQLFTest(glmfit = efit, coef = coef)
  } else if (!is.null(contrast)) {
    if (missing(contrast_matrix)) stop("contast_matrix must be specified when contrast is used")
    qlfTest <- edgeR::glmQLFTest(glmfit = efit, contrast = contrast_matrix[, contrast])
  } else {
    stop("Must specify either coef or contrast")
  }
  out <- edgeR::topTags(qlfTest, n = Inf, p.value = 1)$table
  data.table::setDT(out, keep.rownames = TRUE)
  data.table::setnames(out, "rn", id_name)
  out[]
}
