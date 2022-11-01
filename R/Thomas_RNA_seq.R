ggplot_volcano <- function(df, type="edgeR") {
  if(type=="edgeR") ggplot_volcano_edgeR(df)
  else if(type=="voom")~ggplot_volcano_voom(df)
}

ggplot_volcano_edgeR <- function(df){
  allLogFc <- df[["logFC"]]
  minLogFc <- min(allLogFc)
  maxLogFc <- max(allLogFc)
  minPval <- max(-log10(df[["PValue"]]))
  maxPval <- 0
  ggplot(df, aes(x = logFC, 
                 y = -log10(PValue), 
                 colour = FDR < 0.05)) + 
    geom_point(size=0.5, alpha=0.2) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                       name="Significance", labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR ≥ 0.05"), 
                       breaks = c("TRUE", "FALSE")) +
    ylab(expression(paste(-log[10], "(P-value)"))) +
    xlab(expression(paste(log[2], "(Fold Change)"))) +
    scale_x_continuous(limits = c(minLogFc, maxLogFc)) +
    scale_y_continuous(limits = c(maxPval, minPval)) +
    theme_bw()
}

ggplot_volcano_voom <- function(df){
  allLogFc <- df[["logFC"]]
  minLogFc <- min(allLogFc)
  maxLogFc <- max(allLogFc)
  minPval <- max(-log10(df[["P.Value"]]))
  maxPval <- 0
  p <- ggplot(df, aes(x = logFC, 
                      y = -log10(P.Value), 
                      colour = adj.P.Val < 0.05)) + 
    geom_point(size=0.5, alpha=0.2) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                       name="Significance", 
                       labels = c("TRUE" = "adj.P.Val < 0.05", "FALSE" = "adj.P.Val ≥ 0.05"), 
                       breaks = c("TRUE", "FALSE")) +
    ylab(expression(paste(-log[10], "(P-value)"))) +
    xlab(expression(paste(log[2], "(Fold Change)"))) +
    scale_x_continuous(limits = c(minLogFc, maxLogFc)) +
    scale_y_continuous(limits = c(maxPval, minPval)) +
    theme_bw()
  return(p)
}

#' Return a nicely formatted DT::datatable object of the DE results
#'
#' @param df The dataframe with the DE results
#' @param type Whether the results is edgeR or voom
#'
#' @return A DT::datatable
#' @export
#'
#' @examples
#' get_DE_datatable()
get_DE_datatable <- function(df, type="edgeR"){
  if(type=="edgeR") get_DE_datatable_edgeR(df)
  else if(type=="voom") get_DE_datatable_voom(df)
}

get_DE_datatable_edgeR <- function(df) {
  df %>%
    dplyr::transmute(across(any_of(c("SYMBOL", "ENSEMBL"))),
              description,
              logFC=signif(logFC,2),
              PValue=signif(PValue,2),
              FDR=signif(FDR,2)) |>
    DT::datatable(extensions = 'Buttons', 
                  filter="top", 
                  rownames = FALSE, 
                  options = list(dom = 'Bfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'colvis'),
                                 pageLength = 10))
}

get_DE_datatable_voom <- function(df) {
  df %>% 
    dplyr::transmute(Gene = external_gene_name,
              ENSEMBL,
              description,
              logFC=signif(logFC,2),
              adj.P.Val=signif(adj.P.Val,2),
              P.Value=signif(P.Value,2)) %>% 
    DT::datatable(extensions = 'Buttons', 
                  filter = "top", 
                  rownames = FALSE, 
                  options = list(
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel', 'colvis'),
                    pageLength = 10))
}

#' Return a nicely formatted DT::datatable of the GO-results
#'
#' @param df The GO-results
#'
#' @return A DT::datatable
#' @export
#'
#' @examples
#' get_GO_datatable()
get_GO_datatable <- function(df) {
  df %>% 
    dplyr::transmute(ID, 
              Direction=as.factor(Direction), 
              TERM, 
              `#genes`=NGenes, 
              PValue=signif(PValue,2),
              FDR=signif(FDR, 2)) %>% 
    DT::datatable(extensions = 'Buttons', 
                  filter="top", 
                  rownames = FALSE, 
                  options = list(dom = 'Bfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'colvis'),
                                 pageLength = 10))
}

#' Prints the HTML code to write in the HTML in an asis-chunk, containing the tables of resultss
#'
#' @param df The DE dataframe
#' @param type The type of dataframe
#'
#' @return Raw HTML code
#' @export
#'
#' @examples
#' print_DE_results
print_DE_results <- function(test_list, contrast, type="edgeR") {
  cat("##### Significant genes (FDR<0.05) \n")
  cat("\n")
  test_list[[contrast]] %>% 
    dplyr::filter(FDR<0.05) %>% 
    get_DE_datatable(type = type) %>% 
    knitr:::knit_print() %>% 
    cat()
  cat("\n")
  cat("\n")
  
  
  cat("##### All genes \n")
  cat("\n")
  test_list[[contrast]] %>% 
    get_DE_datatable(type=type) %>%
    knitr:::knit_print() %>% 
    cat()
  cat("\n")
  cat("\n")
  
  cat("##### Volcano plot \n") # Volcano plot
  cat("\n")
  test_list[[contrast]] %>%
    ggplot_volcano(type=type) %>% 
    print()
  cat("\n")
  cat("\n")
}

print_GO_results <- function(ontology_test_list, contrast) {
  cat("##### Significant GO-terms (FDR<0.5) \n")
  cat("\n")
  
  ontology_test_list[[contrast]] %>% 
    filter(FDR<0.05) %>%
    get_GO_datatable() %>% 
    knitr:::knit_print() %>% 
    cat()
  cat("\n")
  cat("\n")
  
  
  cat("##### All GO-terms \n")
  cat("\n")
  ontology_test_list[[contrast]] %>%
    get_GO_datatable() %>%
    knitr:::knit_print() %>% 
    cat()
  cat("\n")
  cat("\n")
}

#' Render one markdown report for a single tissue
#'
#' @param markdown_path The path of the markdown file
#' @param tissue The vector to filter the file by
#' @param exclude The samples to exclude when testing
#'
#' @return Nothing
#' @export
#'
#' @examples
#' render_markdown_file_single_tissue()
render_markdown_file_single_tissue <- function(markdown_path, tissue, exclude) {
  # assuming the output format of input.Rmd is PDF
  rmarkdown::render(
    markdown_path,
    output_file = paste(
      format(Sys.time(),'%Y_%m_%d'),
      tissue,
      sep="_"), 
    # envir=new.env(),
    params = list(tissue = tissue,
                  exclude=exclude),
    envir = parent.frame())
}
#' Render markdown reports for multiple tissues
#'
#' @param markdown_path The path of the template markdown file
#' @param tissue The vector to filter the file by
#' @param exclude The samples to exclude when testing
#'
#' @return Nothing
#' @export
#'
#' @examples
#' render_markdown_file_single_tissue()
render_tissues <- function(markdown_path, tissue_exclusion_vec) {
  for (i in seq_along(tissue_exclusion_vec)) {
    render_markdown_file_single_tissue(markdown_path, 
                                       tissue=names(tissue_exclusion_vec[i]), 
                                       exclude=tissue_exclusion_vec[i])
  }
}

get_enrichment_terms_both <- function (org_db, gene_ids, 
                                       gene_id_type_key_tpe="ENSEMBL", 
                                       min_genes = 5, 
                                       max_genes = 500, 
                                       cache_path) {
  
  if (!(gene_id_type_key_tpe == "ENSEMBL" | gene_id_type_key_tpe == "SYMBOL")) 
    stop("gene_id_type_key_tpe must be either ENSEMBL or SYMBOL")
  
  if (is.character(org_db)) 
    org_db <- get(org_db)
  species_id <- BiocGenerics::species(org_db)
  go_data <- get_go_terms_both(org_db = org_db, 
                               gene_ids = gene_ids, 
                               gene_id_type_key_tpe)
  if (missing(cache_path)) {
    use_cache <- FALSE
    cache_path <- NULL
  }
  reactome_data <- get_reactome_terms_both(org_name = species_id, 
                                           gene_ids = gene_ids, use_cache = use_cache, cache_path = cache_path,
                                           gene_id_type_key_tpe)
  if (nrow(reactome_data) > 0) {
    data.table::set(reactome_data, j = "Species", value = NULL)
    go_data[["Reactome"]] <- reactome_data
  }
  format_enrichment <- function(x) {
    index <- split(x[[gene_id_type_key_tpe]], x[["ID"]])
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

get_go_terms_both <- function (org_db, gene_ids, gene_id_type_key_tpe) 
{
  all_go <- AnnotationDbi::select(org_db, keys = gene_ids, 
                                  columns = "GOALL", keytype = gene_id_type_key_tpe)
  
  data.table::setDT(all_go, key = gene_id_type_key_tpe)
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



get_reactome_terms_both <- function (org_name, gene_ids, use_cache = FALSE, cache_path = NULL, gene_id_type_key_tpe) 
{
  if (use_cache && file.exists(cache_path)) {
    terms <- data.table::fread(cache_path)
  }
  else {
    terms <- data.table::fread("https://reactome.org/download/current/Ensembl2Reactome.txt", 
                               header = FALSE)
    if (!is.null(cache_path)) {
      dir.create(dirname(cache_path), recursive = TRUE, 
                 showWarnings = FALSE)
      data.table::fwrite(terms, cache_path)
    }
  }
  if (use_cache && is.null(cache_path)) 
    warning("use_cache is set but no cache_path is given. No cache is used.")
  if (!missing(org_name)) {
    terms <- terms[org_name, on = "V6", nomatch = FALSE]
  }
  if (!missing(gene_ids)) {
    terms <- terms[gene_ids, on = "V1", nomatch = FALSE]
  }
  . <- V1 <- V2 <- V4 <- V6 <- NULL
  terms <- terms[, .(V1, V2, V4, V6)]
  data.table::setnames(terms, c(gene_id_type_key_tpe, "ID", "TERM", "Species"))
  terms[]
}

#' Plot top 50 DE genes
#'
#' @param y DGE list object
#' @param gene_test_list Output from EdgeR tester
#' @param Character vector containing the variable to order by as present in the y$samples dataframe.
#'
#' @return Plots a heatmap
#' @export
#'
#' @examples
#' plot_top_50_DE_genes_heatmap
plot_top_50_DE_genes_heatmap <- function(y, gene_test_list, order_by_vec) {
  annotation_col_df_DE_heatmap <-  dplyr::select(y$samples, all_of(meta_variables)) %>% 
    magrittr::set_rownames(y$samples$Sample.ID) 
  
  for (i in seq_along(gene_test_list)) {
    top_50_genes <- gene_test_list[[i]] %>% 
      dplyr::filter(FDR<0.05) %>% 
      arrange(FDR) %>% 
      head(50) %>% 
      pull(Gene)
    
    sorted_row_names <- annotation_col_df_DE_heatmap %>% 
      dplyr::arrange(.data[[order_by_vec[i]]]) %>% 
      rownames()
    
    cpm_matrix <- y[c(top_50_genes), sorted_row_names] %>% 
      cpm(log=T)
    print(annotation_col_df_DE_heatmap)
    print(sorted_row_names)
    print(cpm_matrix)
    
    cpm_matrix %>% 
      pheatmap::pheatmap(annotation_col = annotation_col_df_DE_heatmap,
                         display_numbers=F, 
                         cluster_cols = F, 
                         cluster_rows = F,
                         fontsize_row=6 )
  }
}

n_sig_genes_pr_contrast <- function(test_df) {
  
  purrr::map(test_df,
             ~dplyr::filter(., across(any_of(c("FDR", "adj.P.Val")), ~.x < 0.05)) %>% nrow()) %>% 
    unlist()
}