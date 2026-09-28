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

  # the GENCODE v48 GTF is no longer in tximportData (>= 1.41.1), so build
  # a dummy GRanges for the annotated transcripts using the `|`-delimited
  # metadata in the transcript names (IDs and lengths are real, locations are not)
  library(GenomicRanges)
  tnames <- read.delim(files[1])$tname
  tnames <- tnames[grepl("|", tnames, fixed=TRUE)]
  info <- do.call(rbind, strsplit(tnames, "|", fixed=TRUE))
  annotated_gr <- GRanges(
    seqnames = "chrUn",
    ranges = IRanges(start=1, width=as.integer(info[,7])),
    tx_name = info[,1],
    gene_id = info[,2],
    gene_name = info[,6],
    tx_type = info[,8]
  )
  names(annotated_gr) <- annotated_gr$tx_name
  makeLinkedTxpData(
    digest = "6fc626c828b7a342ab0c6ff753055761989bf0e2306370e8766fedf45ad3adb3",
    digestType = "sha256",
    indexName = "gencode.v48",
    txpData = annotated_gr,
    source = "LocalGENCODE",
    organism = "Homo sapiens",
    release = "48",
    genome = "GRCh38"
  )

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
  novel_gr <- as(novel, "GRanges")
  names(novel_gr) <- novel$tx_name
  seqlevels(novel_gr) <- c(seqlevels(novel_gr), seqlevels(annotated_gr))

  # importData for mixed references: first step returns an un-ranged SE
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
  expect_equal(nrow(se_update), nrow(se0))
  expect_equal(sum(mcols(se_update)$index == "annotated", na.rm=TRUE), nrow(info))

  # 22 chr x 500 txps per chrom = 11000 novel txps
  # (these may have metadata if linkedTxpData below was run in a previous session)
  not_in_annotated <- rownames(se_update)[!mcols(se_update)$index %in% "annotated"]
  expect_equal(length(not_in_annotated), 11000L)
  expect_true(all(grepl("novel", not_in_annotated)))

  # `prefer` is respected: linkedTxpData first, even if a linkedTxome
  # for this digest exists in the BFC from a previous session
  prefer <- c("txpdata","txome","precomputed")
  expect_true(inspectDigests(se_mix, prefer=prefer)$linkedTxpData[1])
  se_update <- updateMetadata(se_mix, prefer=prefer)
  # `gene_name` is only in the linkedTxpData, not the TxDb
  expect_equal(mcols(se_update)["ENST00000832824.1", "gene_name"], "DDX11L16")
  # character columns are NA (not "NA") where there is no metadata
  expect_equal(sum(is.na(mcols(se_update)$gene_name)), 11000L)
  expect_error(updateMetadata(se_mix, prefer="foo"))

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
    digestType = "sha256",
    indexName = "my_novel_txps",
    txpData = novel_gr,
    source = "novel", organism="Homo sapiens", 
    release="v1", genome="GRCh38"
  )

  inspectDigests(se_mix)

  inspectDigests(se_mix, count=TRUE)

})

test_that("mergeTxpDataIntoRowData() works with list-types from S4Vectors", {

  # mergeTxpDataIntoRowData() adds missing columns and fills in data for updating rowData
  # initializing list-type vectors from S4Vectors with NA is tested here

  library(S4Vectors)
  rowdata <- DataFrame(test=1:4, row.names=letters[1:4])
  txpDataToAdd <- DataFrame(
    foo=11:12,
    row.names=letters[2:3]
  )
  out <- mergeTxpDataIntoRowData(
    rowdata, txpDataToAdd,
    matches=c("b","c"), indexName="foobar"
  )
  out
  expect_equal(out$foo, c(NA,11,12,NA))

  # again, with list-like vector, `charlist`
  txpDataToAdd2 <- DataFrame(
    charlist = CharacterList("a","c","d"),
    row.names=letters[c(1,3,4)]
  )
  out2 <- mergeTxpDataIntoRowData(
    out, txpDataToAdd2,
    matches=c("a","c","d"), indexName="foobar"
  )
  expect_equal(out2$charlist, CharacterList("a",NA,"c","d"))

  # again, with character and factor vectors (should be NA, not "NA")
  txpDataToAdd3 <- DataFrame(
    char = c("x","y"),
    fac = factor(c("x","y")),
    row.names=letters[2:3]
  )
  out3 <- mergeTxpDataIntoRowData(
    out2, txpDataToAdd3,
    matches=c("b","c"), indexName="foobar"
  )
  expect_identical(out3$char, c(NA,"x","y",NA))
  expect_identical(out3$fac, factor(c(NA,"x","y",NA)))

})
