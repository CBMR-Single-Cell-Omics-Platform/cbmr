#' Get ENSEMBL to Reactome data
#'
#' @param org_name optional string, organism name to subset. Common options are "Mus musculus" and "Homo sapiens"
#' @param ensembl_ids optional character vector. ENSEMBL genes to subset
#' @param cache_path optional path to cache
#'
#' @details Load ENSEMBL to Reactome mappings. If org_name is specified only terms
#' with that species is returned. If ensembl_ids are specified only terms relating to
#' these genes are returned. 
#' 
#' If cache_path is set, the specified file 
#' is loaded if it exists. If it does not exist it is downloaded and the saved to that 
#' location. 
#'
#' @return data.table with ENSEMBL to Reactome mapping
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' get_reactome_terms("Homo sapiens")
#' }
get_reactome_terms <- function (org_name, gene_ids, cache_path = NULL, gene_id_key_type) 
{
  if (!is.null(cache_path) && file.exists(cache_path)) {
    terms <- data.table::fread(cache_path)
  } else {
    terms <- data.table::fread("https://reactome.org/download/current/Ensembl2Reactome.txt", 
                               header = FALSE)
    if (!is.null(cache_path)) {
      dir.create(dirname(cache_path), recursive = TRUE, 
                 showWarnings = FALSE)
      data.table::fwrite(terms, cache_path)
    }
  }
  
  if (!missing(org_name)) {
    terms <- terms[org_name, on = "V6", nomatch = FALSE]
  }
  if (!missing(gene_ids)) {
    terms <- terms[gene_ids, on = "V1", nomatch = FALSE]
  }
  . <- V1 <- V2 <- V4 <- V6 <- NULL
  terms <- terms[, .(V1, V2, V4, V6)]
  data.table::setnames(terms, c(gene_id_key_type, "ID", "TERM", "Species"))
  terms[]
}

#' Get ENSEMBL to GO data
#'
#' @param org_db AnnotationDb class object matching the species analyzed
#' @param ensembl_ids optional character vector. ENSEMBL genes to subset
#' 
#' @details Prepare ENSEMBL to GO mappings. If ensembl_ids are specified only 
#' terms relating to these genes are returned. 
#' 
#' @return list of data.table with ENSEMBL to GO mappings
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' library(org.Mm.eg.db)
#' get_go_terms(org.Mm.eg.db)
#' }
get_go_terms <- function (org_db, gene_ids, gene_id_key_type) 
{
  all_go <- AnnotationDbi::select(org_db, keys = gene_ids, 
                                  columns = "GOALL", keytype = gene_id_key_type)
  
  data.table::setDT(all_go, key = gene_id_key_type)
  if (!missing(gene_ids)) {
    all_go <- all_go[gene_ids, ]
  }
  GOALL <- NULL
  all_go <- all_go[!is.na(GOALL)]
  data.table::setkeyv(all_go, "ONTOLOGYALL")
  data.table::set(all_go, j = "EVIDENCEALL", value = NULL)
  term_go <- AnnotationDbi::select(GO.db::GO.db, keys = AnnotationDbi::keys(GO.db::GO.db), 
                                   columns = c("TERM", "ONTOLOGY"))
  data.table::setDT(term_go)
  all_go_terms <- data.table::merge.data.table(all_go, term_go, 
                                               by.x = c("GOALL", "ONTOLOGYALL"), by.y = c("GOID", "ONTOLOGY"))
  data.table::setnames(all_go_terms, "GOALL", "ID")
  out <- split(all_go_terms, all_go_terms[["ONTOLOGYALL"]])
  out <- lapply(out, data.table::set, j = "ONTOLOGYALL", value = NULL)
  lapply(out, `[`)
}

#' Helper function to subset gene ontology lists
#'
#' @param low minimum number of genes
#' @param high maximum number of genes 
#'
#' @return function to filter lists of genes
make_size_filter <- function(low, high) {
  function(y) {
    x <- length(y)
    (x >= low) & (x <= high)
  }
}

#' Get terms for gene ontology enrichment analysis.
#'
#' @param org_db AnnotationDb class object matching the species analyzed
#' @param ensembl_ids optional character vector. ENSEMBL genes to subset
#' @param min_genes integer, minimum number of genes to include term
#' @param max_genes integer, maximum number of genes to include term
#' @param cache_path character, optional path to cached Reactome file
#'
#' @return list of data formatted to work with the camera and fry wrappers in this package.
#' @export
#'
#' @examples
#' \dontrun{
#' library(org.Mm.eg.db)
#' #rna_counts is an example dataset included in the package
#' get_enrichment_terms(org.Mm.eg.db, rna_counts$Geneid)
#' }
get_enrichment_terms <- function (org_db, gene_ids, 
                                       gene_id_key_type="ENSEMBL", 
                                       min_genes = 5, 
                                       max_genes = 500, 
                                       cache_path) {
  
  if (!(gene_id_key_type == "ENSEMBL" | gene_id_key_type == "SYMBOL")) 
    stop("gene_id_key_type must be either ENSEMBL or SYMBOL")
  
  if (is.character(org_db)) 
    org_db <- get(org_db)
  species_id <- BiocGenerics::species(org_db)
  go_data <- get_go_terms(org_db = org_db, 
                               gene_ids = gene_ids, 
                               gene_id_key_type)
  if (missing(cache_path)) {
    cache_path <- NULL
  }
  reactome_data <- get_reactome_terms(org_name = species_id, 
                                      gene_ids = gene_ids, 
                                      cache_path = cache_path,
                                      gene_id_key_type)
  if (nrow(reactome_data) > 0) {
    data.table::set(reactome_data, j = "Species", value = NULL)
    go_data[["Reactome"]] <- reactome_data
  }
  format_enrichment <- function(x) {
    index <- split(x[[gene_id_key_type]], x[["ID"]])
    index <- lapply(index, unique)
    index <- Filter(cbmr:::make_size_filter(min_genes, max_genes), 
                    index)
    . <- ID <- TERM <- NULL
    annotation <- unique(x[, .(ID, TERM)])
    annotation <- annotation[names(index), on = "ID"]
    list(index = index, annotation = annotation[])
  }
  lapply(go_data, format_enrichment)
}

#' Shallow wrapper to reformat fry/camera output
#'
#' @param contrast name of contrast to be tested or numeric vector of same 
#' length as the number of columns of design. See edgeR::camera.
#' @param index list of index vectors. See edgeR::camera.
#' @param y DGEList object.
#' @param fun function to use for testing. Typically limma::camera or limma::fry
#' @param ... other parameters passed to fun.
#'
#' @return data.table with enrichment test results.
#' @import data.table
enrich_test_wrapper <- function(contrast, index, y, fun, ...)
{
  out <- fun(y, index = index, design = y$design, contrast = contrast, ...)
  data.table::setDT(out, keep.rownames = TRUE)
  data.table::setnames(out, "rn", "ID")
  PValue <- NULL
  out <- out[!is.na(PValue), ]
  out <- out[order(PValue)]
  out[]
}

#' Run ontology tests and annotate results
#'
#' @param terms list with two elements: 'index' containing an index vector as 
#' described in limma::camera and 'annotation' containing information on the 
#' ontologies tested. Can be generated with get_enrichment_terms
#' @param y DGEList object to be tested
#' @param contrast_matrix matrix with contrasts to be tested. Both contrast_matrix 
#' and coefs can be specified simultaneously.
#' @param coefs coefficeints of y$design to be tested Both contrast_matrix 
#' and coefs can be specified simultaneously.
#' @param fun function to use for testing. Typically limma::camera or limma::fry
#' @param ... other parameters passed to fun.
#'
#' @return list of enrichment results
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' library(org.Mm.eg.db)
#' y <- prepare_rna_dgelist(rna_counts, example_metadata, sample_col = "Sample.ID")
#' design <- model.matrix(~ 0 + Genotype + Batch, data = y$samples)
#' y <- edgeR::estimateDisp(y, design = design)
#' rownames(y) <- y$genes$Geneid
#' contrast_matrix <- limma::makeContrasts(Genotype = GenotypeKO - GenotypeWT, 
#'                                         levels = y$design)
#' terms <- get_enrichment_terms(org.Mm.eg.db, y$genes$Geneid)
#' run_ontology_tests(terms$BP, y, contrast_matrix, c("Batch2", "Batch3"), limma::camera)
#' }
run_ontology_tests <- function(terms, y, contrast_matrix = NULL, coefs = NULL, fun, ...) {
  to_test <- list()
  if (!is.null(contrast_matrix)) {
    contrasts <- split(contrast_matrix, col(contrast_matrix))
    names(contrasts) <- colnames(contrast_matrix)
    to_test <- c(to_test, contrasts)
  }
  if (!is.null(coefs)) {
    names(coefs) <- coefs
    to_test <- c(to_test, as.list(coefs))
  }
  if (length(to_test) == 0) stop("contrast_matrix and/or coefs must be specified")
  
  frac_present <- mean(rownames(y) %in% unique(unlist(terms$index)))
  if (frac_present == 0) stop("Rownames of y not present in terms")
  if (frac_present < 1) message(paste0(round(frac_present*100, digits = 2), "% of genes have an annotation"))
  
  enrichments <- lapply(to_test, enrich_test_wrapper, terms$index, y, fun, ...)
  out <- lapply(enrichments, data.table::merge.data.table, y = terms$annotation, 
                by = "ID", nomatch = FALSE)
  PValue <- NULL
  out <- lapply(out, `[`, order(PValue, decreasing = FALSE))
  lapply(out, `[`)
}
