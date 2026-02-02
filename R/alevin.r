# breaking out a lot of alevin specific code that was in main tximeta() body
# Sep 2025 -- this could be cleaned up more...
# mostly this is being preserved as legacy code for importing alevin as SE/SCE
tximetaAlevin <- function(
  coldata,
  type = "alevin",
  txOut,
  skipMeta,
  skipSeqinfo,
  useHub,
  markDuplicateTxps,
  cleanDuplicateTxps,
  customMetaInfo,
  skipFtp,
  ...
) {
  message(paste("importing", type, "quantification files"))

  files <- as.character(coldata$files)

  if (length(files) > 1) {
    stop("alevin import currently only supports a single experiment")
  }

  # remove the files column from colData
  coldata <- subset(coldata, select = -files)

  # tximeta metadata
  metadata <- makeMetadata(type)

  if (skipMeta) {
    txi <- tximport(files, type = type, txOut = txOut, ...)
    metadata$countsFromAbundance <- txi$countsFromAbundance
    coldata <- data.frame(row.names = colnames(txi[["counts"]]))
    se <- makeUnrangedSE(txi, coldata, metadata)
    return(se)
  } else {
    if (!txOut) {
      stop(
        "tximeta is designed to have transcript-level output for Salmon and piscem.
  set txOut=TRUE and use summarizeToGene for gene-level summarization"
      )
    }
  }

  metaInfo <- list(getMetaInfo(
    dirname(files),
    type = "salmon",
    customMetaInfo = customMetaInfo
  ))

  # Salmon's SHA-256 hash of the index is called "index_seq_hash" in the meta_info.json file
  indexSeqHash <- metaInfo[[1]]$index_seq_hash # first sample

  # reshape
  metaInfo <- reshapeMetaInfo(metaInfo, hashType="salmon")
  # add to metadata list
  metadata$quantInfo <- metaInfo

  # try to import files early, so we don't waste user time
  # with metadata magic before a tximport error
  txi <- tximport(files, type = type, txOut = TRUE, ...)
  metadata$countsFromAbundance <- txi$countsFromAbundance

  # try and find a matching txome
  txomeInfo <- getTxomeInfo(digest = indexSeqHash, prefer=c("txome","precomputed"))
  if (is.null(txomeInfo)) {
    message(
      "couldn't find matching transcriptome, returning non-ranged SummarizedExperiment"
    )
    coldata <- data.frame(row.names = colnames(txi[["counts"]]))
    se <- makeUnrangedSE(txi, coldata, metadata)
    return(se)
  }

  # build or load a TxDb from the gtf
  txdb <- getTxDb(txomeInfo, useHub = useHub, skipFtp = skipFtp)

  # build or load transcript ranges (alevin gets gene ranges instead)
  # alevin gets gene ranges instead
  message("generating gene ranges")
  # here gene ranges are named 'txps' for compatibility with code below...
  txps <- getRanges(txdb = txdb, txomeInfo = txomeInfo, type = "gene")
  metadata$level <- "gene"

  # package up the assays
  if ("variance" %in% names(txi)) {
    if ("infReps" %in% names(txi)) {
      assays <- c(txi[c("counts", "variance")], txi$infReps)
      names(assays) <- c(
        "counts",
        "variance",
        paste0("infRep", seq_along(txi$infReps))
      )
    } else {
      assays <- txi[c("counts", "variance")]
    }
  } else {
    assays <- txi["counts"]
  }
  # add mean information as well if it exists in the list
  if ("mean" %in% names(txi)) {
    assays <- c(assays, txi["mean"])
  }
  # add tier information as well if it exists in the list
  if ("tier" %in% names(txi)) {
    assays <- c(assays, txi["tier"])
  }
  coldata <- data.frame(row.names = colnames(assays[["counts"]]))

  # Ensembl FASTA has txp version numbers,
  # but in the Ensembl GTF it is not in the txname,
  # so here we have to remove the version number to build the SummarizedExperiment
  if (txomeInfo$source %in% c("Ensembl", "LocalEnsembl")) {
    txId <- sub("\\..*", "", rownames(assays[["counts"]]))
    for (nm in names(assays)) {
      rownames(assays[[nm]]) <- txId
    }
  }

  # the following functions modifies assays and txps to clean/mark duplicate txps 
  # (this occurs when salmon collapses identical transcripts during indexing)
  if (cleanDuplicateTxps) {
    dup.output.list <- duplicateTxpsClean(
      assays, txps, txomeInfo,
      markDuplicateTxps, cleanDuplicateTxps
    )
    assays <- dup.output.list$assays
    txps <- dup.output.list$txps
  }
  assays <- stripAllCharsAfterBar(assays)
  assays <- checkAssays2Txps(assays, txps)
  txps <- txps[rownames(assays[["counts"]])]
  if (markDuplicateTxps) {
    dup.output.list <- duplicateTxpsMark(
      assays, txps, txomeInfo,
      markDuplicateTxps, cleanDuplicateTxps
    )
    assays <- dup.output.list$assays
    txps <- dup.output.list$txps
  }

  # Ensembl already has nice seqinfo attached...
  # if GENCODE, and not from AHub (which have seqinfo)
  missingSeqinfo <- any(is.na(seqlengths(txps)))
  if (txomeInfo$source == "GENCODE" & !skipSeqinfo & missingSeqinfo) {
    message("fetching genome info for GENCODE")
    ucsc.genome <- genome2UCSC(txomeInfo$genome)
    try(seqinfo(txps) <- Seqinfo(genome = ucsc.genome)[seqlevels(txps)])
  } else if (txomeInfo$source == "RefSeq" & !skipSeqinfo & missingSeqinfo) {
    # if RefSeq...
    message("fetching genome info for RefSeq")
    refseq.genome <- gtf2RefSeq(txomeInfo$gtf, txomeInfo$genome)
    stopifnot(all(seqlevels(txps) %in% seqnames(refseq.genome)))
    try(seqinfo(txps) <- refseq.genome[seqlevels(txps)])
  }

  # add more metadata
  txdbInfo <- metadata(txdb)$value
  names(txdbInfo) <- metadata(txdb)$name
  metadata$txomeInfo <- txomeInfo
  metadata$txdbInfo <- txdbInfo

  se <- SummarizedExperiment(
    assays = assays,
    rowRanges = txps,
    colData = coldata,
    metadata = metadata
  )
  se
}
