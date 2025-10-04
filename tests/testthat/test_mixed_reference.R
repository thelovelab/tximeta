context("mixed reference")
test_that("importing oarfish with mixed reference works as expected", {

  dir <- system.file("extdata/oarfish", package="tximportData")
  names <- paste0("rep", 2:4)
  files <- file.path(dir, paste0("sgnex_h9_", names, ".quant.gz"))
  coldata <- data.frame(files, names)

  # setting: user has run oarfish with, e.g. 
  # --annotated gencode.v48.transcripts.fa.gz 
  # --novel novel.fa.gz

  # skipMeta: we get back the quantification counts but no metadata
  se0 <- tximeta(coldata, type="oarfish", skipMeta=TRUE)

  gtf_dir <- system.file("extdata/gencode", package="tximportData")
  gtf <- file.path(gtf_dir, "gencode.v48.annotation.gtf.gz")
  makeLinkedTxome(
    digest = "6fc626c828b7a342ab0c6ff753055761989bf0e2306370e8766fedf45ad3adb3",
    indexName = "gencode.v48",
    source = "LocalGENCODE",
    organism = "Homo sapiens",
    release = "48",
    genome = "GRCh38",
    fasta = "/path/to/fasta.fa",
    gtf = gtf,
    write = FALSE
  )

  # this prompts them to use importData etc.
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
    end = 1e6 + 1 + 0:499 * 1000 + 1000 - 1,
    strand = "+",
    tx_name = paste0("novel", 1:(22*500)),
    gene_id = paste0("novel_gene", rep(1:(22*10), each=50)),
    type = "protein_coding"
  )
  head(novel)
  library(GenomicRanges)
  novel_gr <- as(novel, "GRanges")
  names(novel_gr) <- novel$tx_name
  seqinfo(novel_gr) <- seqinfo(se)

  # first step just returns an un-ranged SE
  se_mix <- importData(coldata, type="oarfish")
  
  # shows the indices and their digests
  inspectDigests(se_mix)
  # show full digest
  inspectDigests(se_mix, fullDigest=TRUE)
  # this is slower, requires loading the TxDb and ranges...
  inspectDigests(se_mix, count=TRUE)

  # populate what transcript metadata we can find:
  se_update <- updateMetadata(se_mix)
  mcols(se_update)

  # can add ranges, but that requires subsetting to a smaller object 
  # as we can't have a mix of ranges + no-range-data rows
  se_update_w_ranges <- updateMetadata(se_mix, ranges=TRUE)
  mcols(se_update_w_ranges)

  # the user then can add metadata via:
  # linkedTxome() / linkedTxpData() -- they can go do this
  # GRanges or data.frame-like thing
  se_update <- updateMetadata(se_mix, txpData=novel[,-(1:4)])
  mcols(se_update)
  table(mcols(se_update)$index)

  se_update_w_ranges <- updateMetadata(se_mix, txpData=novel_gr, ranges=TRUE)
  mcols(se_update_w_ranges)
  table(mcols(se_update_w_ranges)$index)

  library(BiocFileCache)
  bfc <- BiocFileCache(getBFCLoc())
  bfcinfo(bfc)

  # try out makeLinkedTxpData
  makeLinkedTxpData(
    digest = "43158f2c8e88e3acd77c22aee557625a6f1b6a5038cfc7deb5e64903892d8070",
    indexName = "my_novel_txps",
    txpData = novel_gr,
    source = "novel", organism="Homo sapience", 
    release="v1", genome="GRCh38",
    prefer="txpdata"
  )

})
