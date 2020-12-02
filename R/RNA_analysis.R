#' Read featureCounts output
#'
#' @param x path to a file
#' @param regex optional regular expression for updating sample names. 
#' Default option works with our pipeline. Set to NULL to disable.
#'
#' @return data.table with gene counts and reformatted sample names
merge_featurecounts <- function(files, regex) {
  seqnames <- start <- end <- strand <- NULL
  allCounts <- lapply(files, featurecounts_reader, regex = regex)
  allCounts <- Reduce(f = data.table::merge.data.table, allCounts)
  allCounts[, seqnames:=stringr::str_split_fixed(seqnames, ";", 2)[, 1]]
  allCounts[, start:=as.integer(stringr::str_split_fixed(start, ";", 2)[, 1])]
  allCounts[, end:=as.integer(stringi::stri_reverse(stringr::str_split_fixed(stringi::stri_reverse(end), ";", 2)[, 1]))]
  allCounts[, strand:=stringr::str_split_fixed(strand, ";", 2)[, 1]]
  return(allCounts[])
}

#' Prepare featureCounts output files for RNA-seq analysis
#'
#' @param files list of featureCounts output files
#' @param regex optional regular expression. Default option works with our pipeline
#'
#' @importFrom data.table :=
#' @return data.table of gene counts and formatted gene information. FeatureCounts
#' reports chromosome, start, end and strand for each exon. These are reformatted
#' so start contains the start of the first exon and end the end of the last exon.
#' Only one copy of chromosome and strand are reported. Column names are also 
#' reformatted to match with Bioconductor nomenclature (seqnames instead of Chr, etc.).
#' If a regex is supplied it is used to remove text from the sample names, if not 
#' specified a regex that works with Ali's pipeline is used.
#' @export
#'
#' @examples
#' \dontrun{
#'  featurecount_files <- list.files(path = "path/to", pattern = "files", full.names = TRUE)
#'  count_matrix <- merge_featurecounts(featurecount_files)
#' }
featurecounts_reader <- function(x, regex){
  out <- data.table::fread(x)
  sampleNames <- setdiff(colnames(out), c("Geneid", "Chr", "Start", "End", "Length", "Strand"))
  
  if (missing(regex)){
    sampleNames <- stringr::str_remove(basename(sampleNames), "_S[[:digit:]]+_Aligned.sortedByCoord.out.bam.*")
  } else if (!is.null(regex)){
    sampleNames <- stringr::str_remove_all(sampleNames, regex)
  }
  
  data.table::setnames(out, c("Geneid", "seqnames", "start", "end", "strand", "length", sampleNames))
  data.table::setkeyv(out, c("Geneid", "seqnames", "start", "end", "strand", "length"))
  
  return(out)
}
