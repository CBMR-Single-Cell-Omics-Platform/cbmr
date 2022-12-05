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
  .Deprecated("use catHeader and catHeader_w_tabset instead")

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
  .Deprecated("use catHeader and catHeader_w_tabset instead")
}
