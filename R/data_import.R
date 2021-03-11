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
#' @import data.table
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
#'  featurecount_files <- list.files(path = "path/to", pattern = "txt.gz", full.names = TRUE)
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

#' Prepare a list of cov files for DMR analysis
#'
#' @param files character vector
#' @param BSGenome the BSGenome object for the organism analyzed
#'
#' @return list of data.tables with the same names as files.
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#'  cov_files <- list.files(path = "path/to", pattern = "cov.gz", full.names = TRUE)
#'  count_matrix <- prepare_cov_files(cov_files)
#' }
prepare_cov_files <- function(files, BSGenome = NULL) {
  if (!is.null(BSGenome)) {
    common_chrs <- GenomeInfoDb::standardChromosomes(BSGenome)
    
    count_fn <- function(x){
      cpgs <- Biostrings::matchPattern(pattern = "CG", subject = BSGenome[[x]])
      data.table::data.table(seqnames = x, start = start(cpgs), end = end(cpgs))
    }
    
    all_cpgs <- lapply(common_chrs, count_fn)
    names(all_cpgs) <- common_chrs
    all_cpgs <- data.table::rbindlist(all_cpgs)
    all_cpgs[, id:=seq_len(.N)]
    data.table::setkey(all_cpgs, seqnames, start, end)
  } else {
    all_cpgs <- NULL
  }
  
  lapply(files, covfile_reader, cpgs = all_cpgs)
}

#' Read cov file and optionally aggregate to known CpGs
#'
#' @param file character path to cov file
#' @param cpgs data.table containing the position of known CpGs
#'
#' @return data.table with methylation information, optinally aggregated to known CpGs
#' @import data.table
#' @export
covfile_reader <- function(file, cpgs = NULL) {
  cov_data <- data.table::fread(file, select = c(1:3,5:6))
  data.table::setnames(cov_data, c("seqnames", "start", "end", "Me", "Un"))
  
  if (!is.null(cpgs)) { # Location of known CpGs is supplied
    # Sometimes the reads will use ENSMBL convention (seqnames are e.g. 1, 2, 3)
    # while the BSGenome uses UCSC convention (seqnames are e.g. chr1, chr2, chr3)
    # in that case prepend chr to seqnames
    cpgs_has_chr <- any(stringr::str_detect(cpgs$seqnames, "^chr"))
    seq_has_chr <- any(stringr::str_detect(cov_data$seqnames, "^chr"))
    
    if (cpgs_has_chr & !seq_has_chr){
      cov_data[, seqnames:=paste0("chr", seqnames)]
    }
    
    setkey(cov_data, seqnames, start, end)
    
    ovlp <- foverlaps(cov_data, cpgs)
    
    seqnames <- start <- end <- Me <- Un <- NULL # Silence warnings
    out <- ovlp[, .(seqnames = seqnames[[1]],
             start = start[[1]],
             end = end[[1]],
             Me = sum(Me),
             Un = sum(Un)), by = "id"]
    out[, id:=NULL]
  } else {
    out <- cov_data
  }
  
  out[]
}
