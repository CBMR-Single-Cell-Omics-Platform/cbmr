## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)
rna_link <- "https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm"
rna_homepage <- rvest::read_html(rna_link)
rna_indexes <- rvest::html_table(rna_homepage)[[1]][, c(1, 2, 6)]
data.table::setDT(rna_indexes)
data.table::setnames(rna_indexes, c("index", "i7", "i5"))

rrbs_16 <- data.table::fread("https://www.nugen.com/sites/default/files/Barcodes%3A_Ovation_RRBS_Methyl-Seq_Systems_1-16_2232.txt", skip = 2)
data.table::setnames(rrbs_16, c("index", "i7"))
rrbs_16[, i5:=NA_character_]


rrbs_96 <- data.table::fread("https://www.nugen.com/sites/default/files/Barcodes%3A_Ovation_RRBS_Methyl-Seq_Systems_1-96_2233.txt", skip = 2)
data.table::setnames(rrbs_96, c("index", "del", "i7"))
rrbs_96[, del:=NULL]
rrbs_96[, i5:=NA_character_]


bioo <- readxl::read_excel(here::here("data-raw/Bioo-Scientific-Small-RNA-Barcode-Indices-v1-1-15-1.xlsx"), skip = 2)
data.table::setDT(bioo)
data.table::setnames(bioo, c("index", "i7", "del"))
bioo[, del:=NULL]
bioo <- bioo[1:48, ]
bioo[, i5:=NA_character_]

pki <- readxl::read_excel(here::here("data-raw/PKI-small-RNA-v3-w-UDI-Sequences-v1-30-19.xlsx"), skip = 2)
data.table::setDT(pki)
data.table::setnames(pki, c("index", "i7", "i5"))
pki <- pki[1:192, ]

indexes <- list(
  "Illumina–TruSeq DNA and RNA UD Indexes" = rna_indexes, 
  "NuGEN Ovation RRBS Methyl-Seq System 1-16" = rrbs_16, 
  "NuGEN Ovation RRBS Methyl-Seq System 1-96" = rrbs_96, 
  "Bioo NEXTflex Small RNA Barcodes" = bioo, 
  "PKI NEXTFLEX small RNA-seq v3 with UDI primers" = pki
  )
lapply(indexes, data.table::setkey, index)
usethis::use_data(indexes)
