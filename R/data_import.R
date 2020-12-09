#' Read featureCounts output
#'
#' @param x path to a file
#' @param regex optional regular expression for updating sample names. 
#' Default option works with the fq2count pipeline. Set to NULL to disable.
#'
#' @return data.table with gene counts and reformatted sample names
featurecounts_reader <- function(x, regex){
  out <- data.table::fread(x)
  sampleNames <- setdiff(colnames(out), c("Geneid", "Chr", "Start", "End", "Length", "Strand"))
  
  if (missing(regex)){
    sampleNames <- stringr::str_remove_all(basename(sampleNames), "_S[[:digit:]]+_Aligned.sortedByCoord.out.bam.*|SCOP_[[:digit:]]{4}_")
  } else if (!is.null(regex)){
    sampleNames <- stringr::str_remove_all(sampleNames, regex)
  }
  
  data.table::setnames(out, c("Geneid", "seqnames", "start", "end", "strand", "length", sampleNames))
  data.table::setkeyv(out, c("Geneid", "seqnames", "start", "end", "strand", "length"))
  
  out
}

#' Prepare featureCounts output files for RNA-seq analysis
#'
#' @param files list of featureCounts output files
#' @param regex optional regular expression for updating sample names. 
#' Default option works with the fq2count pipeline. Set to NULL to disable.
#'
#' @importFrom data.table :=
#' @return data.table of gene counts and formatted gene information. FeatureCounts
#' reports chromosome, start, end and strand for each exon. These are reformatted
#' so start contains the start of the first exon and end the end of the last exon.
#' Only one copy of chromosome and strand are reported. Column names are also 
#' reformatted to match with Bioconductor nomenclature (seqnames instead of Chr, etc.).
#' If a regex is supplied it is used to remove text from the sample names, if not 
#' specified a regex that works with the fq2count pipeline is used.
#' @export
#'
#' @examples
#' \dontrun{
#'  featurecount_files <- list.files(path = "path/to", pattern = "files", full.names = TRUE)
#'  count_matrix <- prepare_featurecounts(featurecount_files)
#' }
prepare_featurecounts <- function(files, regex) {
  seqnames <- start <- end <- strand <- NULL
  allCounts <- lapply(files, featurecounts_reader, regex = regex)
  allCounts <- Reduce(f = data.table::merge.data.table, allCounts)
  allCounts[, seqnames:=stringr::str_split_fixed(seqnames, ";", 2)[, 1]]
  allCounts[, start:=as.integer(stringr::str_split_fixed(start, ";", 2)[, 1])]
  allCounts[, end:=as.integer(stringi::stri_reverse(stringr::str_split_fixed(stringi::stri_reverse(end), ";", 2)[, 1]))]
  allCounts[, strand:=stringr::str_split_fixed(strand, ";", 2)[, 1]]
  
  data.table::setkeyv(allCounts, c("Geneid", "seqnames", "start", "end", "strand", "length"))
  
  allCounts[]
}


#' Prepare a DGEList object
#'
#' @param counts data.table with gene counts, if data.table has keys those are used to determine gene information 
#' @param metadata data.table with metadata on the samples
#' @param sample_col column in metadata that corresponds to sample (column) names in counts
#'
#' @return DGEList object with information on both genes and samples
#' @export
#'
#' @examples
#' \dontrun{
#'   prepare_rna_dgelist(counts, metadata, "Sample ID")
#' }
prepare_rna_dgelist <- function(counts, metadata, sample_col = "Sample ID") {
  if (data.table::haskey(counts)) {
    genes <- data.table::key(counts)
    sample_ids <- setdiff(colnames(counts), genes)
  } else {
    genes <- c("Geneid", "seqnames", "start", "end", "strand", "length")
    if (all(genes %in% colnames(counts))) {
      sample_ids <- setdiff(colnames(counts), genes)
    } else {
      stop("counts must either be keyed or contain the columns Geneid, seqnames, start, end, strand & length")
    }
  }
  
  metadata_samples <- metadata[[sample_col]]
  
  if (!setequal(metadata_samples, sample_ids)){
    stop("The same samples are not present in the count matrix and the metadata")
  }
  
  y <- edgeR::DGEList(counts  = counts[, ..metadata_samples], 
                      genes   = counts[, ..genes], 
                      samples = metadata)
  y
}

