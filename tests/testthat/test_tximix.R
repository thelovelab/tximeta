context("tximix")
test_that("tximix works as expected", {

  dir <- system.file("extdata/oarfish", package="tximportData")
  names <- paste0("rep", 2:4)
  files <- file.path(dir, paste0("sgnex_h9_", names, ".quant.gz"))
  coldata <- data.frame(files, names)

  # skipMeta: we get back the quantification counts but no metadata
  se0 <- tximeta(coldata, type="oarfish", skipMeta=TRUE)

  # try to import metadata: it will look only at the annotated digest
  expect_warning({
    se <- tximeta(coldata, type="oarfish")
  }, "the annotation is missing some transcripts")

  not_in_annotated <- rownames(se0)[!rownames(se0) %in% rownames(se)]
  
  # 22 chr x 500 txps per chrom = 11000 novel txps
  expect_equal(sum(grepl("novel",not_in_annotated)), 11000L)

  # rowData(se) # has tx_id, gene_id, tx_name from TxDb also ranges

  # define novel set so we can add metadata
 novel <- data.frame(
    seqnames = paste0("chr", rep(1:22, each=500)),
    start = 1e6 + 1 + 0:499 * 1000,
    width = 1000, strand = "+",
    tx_id = paste0("novel", 1:(22*500)),
    gene_id = paste0("novel_gene", rep(1:(22*10), each=50)),
    type = "protein_coding"
  )
  novel$end <- novel$start + novel$width - 1
  head(novel)
  library(GenomicRanges)
  novel_gr <- as(novel, "GRanges")
  seqinfo(novel_gr) <- seqinfo(se)

  # first step just returns an unranged SE
  se_mix <- tximix(coldata, type="oarfish")
  
  # shows the indices and their digests
  tximixInspectDigests(se_mix)

  # populate what transcript metadata we can find:
  se_update <- tximixUpdateTxpData(se_mix)
  mcols(se_update)

  # maybe then the user wants to add metadata via:
  # linkedTxome -- they can go do this
  # linkedTxpData -- they can go do this
  # GRanges
  # data.frame
  tximixAddTxpData(se_mix, txpData=novel)

})
