#' Title
#'
#' @param file path to project
#' @param ... 
#'
#' @return
#' @export
rna_seq_project <- function(file, ...){
  # ensure path exists
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  here::i_am(path)
  
  folders <- file.path(path, c("data", "data/raw", "R", "output"))
  lapply(folders, dir.create, recursive = TRUE, showWarnings = FALSE)

  # create files
  rmarkdown::draft(file = here::here(paste0("output/", project_name, ".Rmd")), 
                   template = "inst/rmarkdown/templates/bulk-rna-seq", 
                   edit = FALSE)
  writeLines(counts_r, con = here::here(paste0("data/counts.R")))
  writeLines(metadata_r, con = here::here(paste0("data/metadata.R")))
  
  # prepare git (most is taken from the usethis::use_git, 
  # but unfortunately some functions are not exported from usethis)
  gert::git_init(git_repo, bare = TRUE)
  gert::git_init(path)
  
  gert::git_remote_add(git_repo)
  gert::git_fetch(remote = git_repo)
  
  git_ignore_files <- c(
    ".Rproj.user",
    ".Rhistory",
    ".Rdata",
    ".httr-oauth",
    ".DS_Store",
    "*xlsx",
    "*smk",
    "*multiqc*",
    "*csv.gz"
  )
  usethis::write_union(file.path(path, ".gitignore"), git_ignore_files)
  
  if (rstudioapi::hasFun("executeCommand")) {
    rstudioapi::executeCommand("vcsRefresh")
  }
  
  gert::git_add(".")
  
  gert::git_commit(message = "Initial Commit")
  gert::git_push()
}

counts_r <- "library(org.XX.eg.db)
library(data.table)
files <- list.files(here::here('data/raw/featureCounts'),
                    pattern = '.gz',
                    full.names = TRUE
)
counts <- cbmr::prepare_featurecounts(files)
counts[, seqnames := paste0('chr', seqnames)]

conv <- clusterProfiler::bitr(counts$Geneid,
                              fromType = 'ENSEMBL',
                              toType = c('SYMBOL', 'GENENAME'), 
                              OrgDb = org.XX.eg.db
)
data.table::setDT(conv, key = 'ENSEMBL')
conv <- conv[counts$Geneid, , mult = 'first']

counts <- data.table::merge.data.table(
  x = counts,
  y = conv,
  by.x = 'Geneid',
  by.y = 'ENSEMBL'
)
data.table::setcolorder(counts, c(
  'Geneid', 'seqnames', 'start', 'end',
  'strand', 'length', 'SYMBOL', 'GENENAME'
))
data.table::setnames(counts, 'Geneid', 'ENSEMBL')
data.table::setkeyv(counts, cols = c(
  'ENSEMBL', 'seqnames', 'start', 'end',
  'strand', 'length', 'SYMBOL', 'GENENAME'
))
usethis::use_data(counts, overwrite = TRUE)
"

metadata_r <- "metadata <- readxl::read_excel(here::here(
  'data', 'raw',
  'METADATA_FILE.xlsx'
))
idx_filled <- lapply(metadata, function(x) !all(is.na(x)))
metadata <- metadata[, unlist(idx_filled)]
data.table::setDT(metadata)
data.table::setnames(metadata, stringr::str_replace_all(colnames(metadata), 
                                                        '[:blank:]', ' '))
data.table::setnames(
  metadata,
  c(
    'SCOP #', 'Subject Name (User selected ID)', 'Sample ID (Library ID)', 
    'Project name', 'Nucleic acid', 'Species', 'Cell type', 
    'Extraction protocol', 'RIN value', 'Library protocol', 'Library Pool', 
    'Sequencing Lane', 'Run Folder', 'Sequencer', 'Read length', 'Read layout', 
    'Index 1', 'sequence Index 1', 'Index 2', 'sequence Index 2', 'Comments'
  ),
  c(
    'scop_id', 'user_id', 'sample_id', 'project_name',
    'nucleic_acid', 'species', 'cell_type', 'extraction', 'rin',
    'library_protocol', 'library_pool', 'sequencing_lane',
    'run_folder', 'sequencer', 'read_length', 'read_layout',
    'index_1', 'index_1_sequence', 'index_2',
    'index_2_sequence', 'comment'
  )
)

usethis::use_data(metadata, overwrite = TRUE)"