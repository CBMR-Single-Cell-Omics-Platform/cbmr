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