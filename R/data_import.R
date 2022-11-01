#' Prepare metadata data.tale from excel file
#'
#' @param file string, path to excel file
#' @param new_col_names optional named character vector, see details
#' @param uninteresting_cols optional character vector, see details
#' @param use_defult_new_names boolean, use the default name conversions found 
#' in the variable computer_friendly_names?
#' @param use_default_uninteresting boolean, use the default list of 
#' uninteresting columns found in the variable uninteresting_names?
#'
#' @details 
#' new_col_names is a named vector of the format 
#' c("excel col name" = "updated_name"). An example is included in the variable
#' computer_friendly_names.
#' 
#' uninteresting_cols is a vector of columns that will be removed. This removal
#' happens after column names are updated, so if a column is updated with 
#' new_col_names, the updated name needs to be given in uninteresting_cols.
#' A default list is included in the variable uninteresting_names.
#' 
#' @return data.table with formatted metadata
#' @import data.table
#' @export
prepare_metadata <- function(file, new_col_names = NULL,
                             uninteresting_cols = NULL,
                             use_defult_new_names = TRUE,
                             use_default_uninteresting = TRUE) {
  metadata <- readxl::read_excel(file)
  data.table::setDT(metadata)
  metadata <- metadata[, !apply(metadata, 2, function(x) all(is.na(x))), with = FALSE]
  
  conv <- NULL
  if (!is.null(new_col_names)) {
    conv <- c(conv, new_col_names)
  }
  if (use_defult_new_names) {
    conv <- c(conv, computer_friendly_names)
  }
  
  to_remove <- NULL
  if (!is.null(uninteresting_cols)) {
    to_remove <- c(to_remove, uninteresting_cols)
  }
  if (use_default_uninteresting) {
    to_remove <- c(to_remove, uninteresting_names)
  }
  
  if (!is.null(conv)) {
    data.table::setnames(metadata,
             old = names(conv),
             new = conv,
             skip_absent = TRUE)
  }
  
  if (!is.null(to_remove)) {
    to_remove <- setdiff(colnames(metadata), to_remove)
    metadata <- metadata[, .SD, .SDcols = to_remove]
  }
  
  if ("rin" %in% colnames(metadata)) {
    rin <- NULL
    metadata[, rin := as.numeric(rin)]
  }
  
  metadata[]
}


#' Read featureCounts output
#'
#' @param x path to a file
#' @param regex optional regular expression for updating sample names.
#' Default option works with the fq2count pipeline. Set to NULL to disable.
#'
#' @return data.table with gene counts and reformatted sample names
featurecounts_reader <- function(x, regex = NULL){
  out <- data.table::fread(x)
  sampleNames <- setdiff(colnames(out), c("Geneid", "Chr", "Start", "End", "Length", "Strand"))

  if (is.null(regex)){
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
#' @param parallel logical, run in parallel
#' Default option works with the fq2count pipeline. Set to NULL to disable.
#'
#' @import data.table
#' @importFrom foreach %dopar%
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
#' featurecounts_files <- c(
#'   system.file("extdata", "KO_1.txt.gz", package = "cbmr"),
#'   system.file("extdata", "WT_1.txt.gz", package = "cbmr"),
#'   system.file("extdata", "KO_2.txt.gz", package = "cbmr"),
#'   system.file("extdata", "WT_2.txt.gz", package = "cbmr"),
#'   system.file("extdata", "KO_3.txt.gz", package = "cbmr"),
#'   system.file("extdata", "WT_3.txt.gz", package = "cbmr")
#' )
#' rna_counts <- prepare_featurecounts(featurecounts_files)
#' }
prepare_featurecounts <- function(files, regex = NULL, parallel = TRUE) {
  seqnames <- start <- end <- strand <- NULL

  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(cores[1] - 1)
  doParallel::registerDoParallel(cl)

  allCounts <- foreach::foreach(file = files) %dopar% {featurecounts_reader(file, regex = regex)}

  parallel::stopCluster(cl)

  #allCounts <- lapply(files, featurecounts_reader, regex = regex)
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
#'   prepare_rna_dgelist(counts = rna_counts,
#'                       metadata = example_metadata,
#'                       sample_col = "Sample.ID"
#'                       )
#' }
prepare_rna_dgelist <- function(counts, metadata, sample_col = "Sample.ID") {
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
  
  # as.character makes sure samples named 1,2,3... etc select the correct columns
  metadata_samples <- as.character(metadata[[sample_col]])

  if (!setequal(metadata_samples, sample_ids)){
    stop("The same samples are not present in the count matrix and the metadata")
  }

  y <- edgeR::DGEList(counts  = counts[, metadata_samples, with = FALSE],
                      genes   = counts[, genes, with = FALSE],
                      samples = metadata)
  y
}

#' Prepare a list of cov files for DMR analysis
#'
#' @param files character vector
#' @param BSGenome the BSGenome object for the organism analyzed
#' @param min_reads integer, minimum number of reads to include position.
#' Set to 0 to retain all data. Defaults to 8
#' @param verbose Print progress updates
#'
#' @description This helper function reads in Bismark cov(.gz) files. If a BSGenome
#' object is supplied it aggregates the C's on opposite strands to a single CpG
#' measurement, and retains only canonical CpGs. This function is very slow,
#' even on the small example shown below. For a full size experiment runtime
#' is likely 15 - 30 minutes.
#'
#' The function attempt to correctly handle situation where the cov files follow
#' the ENSEMBL scheme of naming chromosomes 1, 2, ... instead of chr1, chr2, ...
#'
#' @return list of data.tables with the same names as files.
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' require(BSgenome.Mmusculus.UCSC.mm10)
#'
#' cov_files <- c(
#'   system.file("extdata", "KO_1.cov.gz", package = "cbmr"),
#'   system.file("extdata", "WT_1.cov.gz", package = "cbmr"),
#'   system.file("extdata", "KO_2.cov.gz", package = "cbmr"),
#'   system.file("extdata", "WT_2.cov.gz", package = "cbmr"),
#'   system.file("extdata", "KO_3.cov.gz", package = "cbmr"),
#'   system.file("extdata", "WT_3.cov.gz", package = "cbmr")
#' )
#'
#' names(cov_files) <- gsub(".cov.gz$", "", basename(cov_files))
#'
#' cov_list <- prepare_cov_files(files = cov_files,
#'                               BSGenome = BSgenome.Mmusculus.UCSC.mm10
#' )
#' }
prepare_cov_files <- function(files, BSGenome = NULL, min_reads = 8, verbose = TRUE) {
  start <- end <- id <- seqnames <- .N <- NULL

  if (!is.null(BSGenome)) {
    if (verbose) message("Preparing CpG map")
    common_chrs <- GenomeInfoDb::standardChromosomes(BSGenome)

    count_fn <- function(x){
      cpgs <- Biostrings::matchPattern(pattern = "CG", subject = BSGenome[[x]])
      data.table::data.table(seqnames = x,
                             start = BiocGenerics::start(cpgs),
                             end = BiocGenerics::end(cpgs))
    }

    all_cpgs <- lapply(common_chrs, count_fn)
    names(all_cpgs) <- common_chrs
    all_cpgs <- data.table::rbindlist(all_cpgs)
    all_cpgs[, id:=seq_len(.N)]
    data.table::setkey(all_cpgs, seqnames, start, end)
  } else {
    all_cpgs <- NULL
  }
  out <- list()
  for (i in seq_along(files)) {
    if (verbose) message("Loading file ", i, " of ", length(files))
    out[[i]] <- covfile_reader(files[[i]], min_reads = min_reads, cpgs = all_cpgs)
  }
  names(out) <- names(files)
  out
}

#' Read cov file and optionally aggregate to known CpGs
#'
#' @param file character path to cov file
#' @param cpgs data.table containing the position of known CpGs
#' @param min_reads minimum number of reads needed to retain CpG
#'
#' @return data.table with methylation information, optionally aggregated to known CpGs
#' @import data.table
#' @export
covfile_reader <- function(file, min_reads = 8, cpgs = NULL) {
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

    seqnames <- start <- end <- Me <- Un <- . <- id <-  NULL # Silence warnings
    out <- ovlp[, .(seqnames = seqnames[[1]],
             start = start[[1]],
             end = end[[1]],
             Me = sum(Me),
             Un = sum(Un)), by = "id"]
    out[, id:=NULL]
  } else {
    out <- cov_data
  }

  out[(Me + Un) > min_reads & !is.na(start), ]
}

#' Remove samples with few sites covered
#'
#' @param cov_data list of data.tables returned by prepare_cov_files
#' @param min_sites minimum number of sites covered
#'
#' @return list of data.tables with min_sites or more covered CpGs
#' @export
remove_shallow_samples <- function(cov_data, min_sites) {
  filter_fn <- function(x) nrow(x) >= min_sites
  Filter(f = filter_fn, x = cov_data)
}

#' Merge list of cov data
#'
#' @param cov_data list of data.tables returned by prepare_cov_files. Elements
#' in cov_data are updated by reference when running this function.
#' Elements must be named. Names are used to construct column names in the
#' output.
#'
#' @return data.table with the rows present in all samples
#' @export
merge_cov <- function(cov_data) {
  if (is.null(names(cov_data))) stop("Elements in cov_data must be named")

  for (i in names(cov_data)) {
    data.table::setnames(cov_data[[i]], c("Me", "Un"), paste0(i, c("_Me", "_Un")))
  }

  merger <- function(x, y) {
    data.table::merge.data.table(
      x = x, y = y, by = c("seqnames", "start", "end")
      )
  }
  Reduce(merger, cov_data)
}

#' Filter methylation data
#'
#' @param x data.table with methylation information.
#' @param remove optional data.table with regions to be removed
#' @param select optional data.table with regions to be selected
#'
#' @details all three arguments must contain columns named seqnames, start and end.
#' If regions are present in both remove and select, they are removed.
#'
#' @return data.table filtered based on the tables in remove and select
#' @export
filter_methylation <- function(x, remove = NULL, select = NULL) {
  if (!is.null(remove)) {
    data.table::setkeyv(remove, c("seqnames", "start", "end"))

    idx_remove <- data.table::foverlaps(x, remove, which = TRUE, nomatch = NULL)
    x <- x[-unique(idx_remove$xid), ]

  }
  if (!is.null(select)) {
    data.table::setkeyv(select, c("seqnames", "start", "end"))

    idx_select <- data.table::foverlaps(x, select, which = TRUE, nomatch = NULL)
    x <- x[unique(idx_select$xid), ]
  }

  x
}

#' Annotate methylation data
#'
#' @param x data.table with methylation information.
#' @param annotation data.table with annotation information
#' @param mult how to handle multi overlaps. See [data.table::foverlaps]
#' for details.
#' @param annotation_col optional, column name of annotation to add to x. If
#' null (the default) the entire table is added.
#' @param keep_annotation optional logical, should the start and end columns of
#' the annotation be retained? If TRUE, these columns are named feature_start and
#' feature_end, or, if annotation_col is specified, that name is used instead of
#' feature.
#'
#' @details x and annotation arguments must contain columns named seqnames,
#' start and end. If the elements in annotation overlap, the overlaps are handled
#' using the mult argument. Defaults to "first" which means only the first feature
#' in annotation is used. Alternatively "all" can be used, but that leads to a
#' CpG being tested multiple times.
#'
#' @return data.table x annotated with the information in annotation
#'
#' @export
annotate_methylation <- function(x, annotation, mult = "first",
                                 annotation_col = NULL, keep_annotation = FALSE) {

  if (!is.null(annotation_col)) {
    if (length(annotation_col) > 1) stop("Only one column can be selected using annotation.")
    if (!(annotation_col %in% colnames(annotation))) stop(paste(annotation_col, "not in column names of annotation."))
    data.table::setkey(annotation, "seqnames", "start", "end")

    annotation <- annotation[, c("seqnames", "start", "end", annotation_col), with = FALSE]
    feature_name <- annotation_col
  } else {
    feature_name <- "feature"
  }

  x <- data.table::foverlaps(x, annotation,
                                    by.x = c("seqnames", "start", "end"),
                                    by.y = c("seqnames", "start", "end"),
                                    type = "any",
                                    mult = mult)
  if(keep_annotation) {
    data.table::setnames(x, c("start", "end"), paste0(feature_name, c("_start", "_end")))
  } else {
    x[, c("start", "end"):=NULL]
  }

  data.table::setnames(x, c("i.start", "i.end"), c("start", "end"))
  x[]
}

#' Wrapper function for preprocessing methylation data
#'
#' @param x data.table with methylation information.
#' @param remove optional, list of files with regions to remove
#' @param select optional, list of files with regions to select
#' @param genes optional, gtf file containing gene information
#' @param annotations optional, list of files containing additional annotation information.
#' If the list is named, only the column with that name will be used. Mixing named
#' and unnamed elements is supported.
#'
#' @details The files in remove, select and annotation must be
#' (optionally gzipped) csv files, or one of the filetypes supported by \link[rtracklayer]{import}.
#' If file is a csv file, the first three columns must be seqnames, start and end.
#'
#' @return
#' @export
preprocess_methylation_data <- function(x, remove = NULL, select = NULL,
                                        genes = NULL, annotations = NULL) {
  importer <- function(file) {
    if (stringr::str_ends(file, "\\.csv(\\.gz)?$")) {
      out <- data.table::fread(file = file)
    } else {
      out <- rtracklayer::import(con = file)
      out <- data.table::as.data.table(out)
    }
    out[]
  }

  collapse <- function(x) {
    x <- lapply(x, function(x) x[, c("seqnames", "start", "end"), with = FALSE])
    out <- do.call(what = "rbind", x)
    out[]
  }

  if (!is.null(remove)) {
    remove <- lapply(remove, importer)
    remove <- collapse(remove)
  }

  if (!is.null(select)) {
    select <- lapply(select, importer)
    select <- collapse(select)
  }

  x <- filter_methylation(x, remove = remove, select = select)

  annot_list <- list()

  if (!is.null(genes)) {
    genes <- rtracklayer::import(genes)

    exons <- genes[genes$type == "exon" & genes$gene_type == "protein_coding"]
    exons <- data.table::as.data.table(exons)
    exons <- exons[, c("seqnames", "start", "end", "gene_id"), with = FALSE]
    data.table::setkeyv(exons, c("seqnames", "start", "end"))
    data.table::setnames(exons, "gene_id", "gene_exon_id")

    promoters <- genes[genes$type == "gene" & genes$gene_type == "protein_coding"]
    promoters <- IRanges::promoters(promoters, upstream = 3000, downstream = 1000)
    promoters <- data.table::as.data.table(promoters)
    promoters <- promoters[, c("seqnames", "start", "end", "gene_id"), with = FALSE]
    data.table::setkeyv(promoters, c("seqnames", "start", "end"))
    data.table::setnames(promoters, "gene_id", "promoter_id")

    annot_list[["gene_exon_id"]] = exons
    annot_list[["promoter_id"]] = promoters
  }

  if (!is.null(annotations)) {
    annot_list <- c(annot_list, lapply(annotations, importer))
  }

  if (length(annot_list) > 0) {
    for (i in seq_len(length(annot_list))) {
      col_name <- names(annot_list)[[i]]
      if (col_name == "") col_name <- NULL
      x <- annotate_methylation(x = x, annotation = annot_list[[i]],
                                annotation_col = col_name)
    }
  }

  x[]
}

#' Prepare a DGEList object for RRBS analysis
#'
#' @param counts data.table with gene counts
#' @param metadata data.table with metadata on the samples
#' @param sample_col column in metadata that corresponds to sample (column) names in counts
#'
#' @return DGEList object with information on both genes and samples
#' @importFrom data.table .SD
#' @export
prepare_rrbs_dgelist <- function(counts, metadata, sample_col = "Sample.ID", aggregate_by = NULL) {
  if (!(sample_col %in% colnames(metadata))) stop(paste(sample_col, "not found in metadata"))

  count_idx <- paste0(rep(metadata[[sample_col]], each = 2), c("_Me", "_Un"))

  if (!all(count_idx %in% colnames(counts))) {
    warning("Not all samples in metadata are present in counts. Subsetting to only
            include samples in count. Make sure this is intended.")
    count_idx <- count_idx[count_idx %in% colnames(counts)]
  }

  if (!is.null(aggregate_by)) {
    counts <- counts[, lapply(.SD, sum), by = aggregate_by, .SDcols = patterns("_Me$|_Un$")]
  }

  genes_idx <- setdiff(colnames(counts), count_idx)

  y <- edgeR::DGEList(counts  = counts[, count_idx, with = FALSE],
                      genes   = counts[, genes_idx, with = FALSE]
                      )
  y
}
