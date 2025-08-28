context("tximix")
test_that("tximix works as expected", {

  dir <- system.file("extdata/oarfish", package="tximportData")
  names <- paste0("rep", 2:4)
  files <- file.path(dir, paste0("sgnex_h9_", names, ".quant.gz"))
  coldata <- data.frame(files, names)
  se0 <- tximeta(coldata, type="oarfish", skipMeta=TRUE)

  se <- tximeta(coldata, type="oarfish")

  not_in_annotated <- rownames(se0)[!rownames(se0) %in% rownames(se)]
  table(grepl("novel",not_in_annotated))

  rowData(se)

  novel <- data.frame(
    seqnames = rep(1:22, each=500),
    start = 1e6 + 1 + 0:499 * 1000,
    width = 1000, strand = "+",
    id = paste0("novel", 1:(22*500))
  )
  novel$end <- novel$start + novel$width - 1
  library(GenomicRanges)
  novel <- as(novel, "GRanges")

  se <- tximix(coldata, type="oarfish")
  
})
