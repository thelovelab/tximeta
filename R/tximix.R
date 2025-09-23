#' Import transcript quantification across sets of transcripts
#'
#' tximix allows separation of annotated and novel transcripts
#' during transcript quantification import and addition of metadata.
#' This function supports quantification in the case of mixed
#' transcriptomes, including additions from de novo assembly,
#' as well as transgenes and spike-ins.
#'
#' @param coldata data.frame with columns \code{files} and \code{names}
#' as in \code{\link{tximeta}}
#' @param type what quantifier was used (see \code{\link{tximport}})
#' @param quiet whether to suppress messages
#' @param ... passed to tximport
#'
#' @return an unranged SummarizedExperiment
#'
#' @export
tximix <- function(coldata, type="oarfish", quiet=FALSE, ...) {
  stopifnot(type == "oarfish")
  
  # tximeta metadata
  metadata <- makeMetadata(type)

  files <- as.character(coldata$files)
  names(files) <- coldata$names
  txi <- tximport(files, type=type, txOut=TRUE, ...)
  metadata$countsFromAbundance <- txi$countsFromAbundance

  # `metaInfo` = list with quantification tool metadata from JSON files
  # that are alongside quantification files in newer tools
  metaInfo <- lapply(
    files,
    getMetaInfo,
    type=type
  )

  # different styles of storing hash value by method
  hashType <- type2hashType(type)
  
  # check the sequence digest (hash) of the transcriptome index with 1st sample
  # readIndexSeqHash() returns a list of functions
  # TODO - Sep 2025: currently this code just looks at the annotated digest,
  # could also look at the novel and confirm it is also consistent
  indexSeqHash <- readIndexSeqHash()[[hashType]](metaInfo[[1]])
  if (length(files) > 1) {
    hashes <- sapply(metaInfo, readIndexSeqHash()[[hashType]])
    if (!all(hashes == indexSeqHash)) {
      stop("the samples do not share the same index, and cannot be imported")
    }
  }
  
  # reshape this list object, invert the JSON hierarchy 
  # and examine consistency of the digest 'index_seq_hash'
  metaInfo <- reshapeMetaInfo(metaInfo, hashType="oarfish")
  
  # add the per-sample metadata from quantification JSON files to the metadata list object
  metadata$quantInfo <- metaInfo

  # assemble list of matrices for outputting an unranged SE
  assays <- txi[c("counts","abundance","length")]

  # GENCODE usually has characters after the ENST... 
  # these disrupt metadata operations (adding ranges or IDs)
  assays <- stripAllCharsAfterBar(assays)

  se <- makeUnrangedSE(assays, coldata, metadata)
  
  if (!quiet)
    message("returning unranged SummarizedExperiment, other tximix functions:\n",
            "-- tximixInspectDigests() to check matching digests\n",
            "-- linkedTxome() / linkedTxpData() to link new digests to GTF / custom metadata\n",
            "-- tximixUpdateTxpData() to update metadata and optionally add ranges"
          )

  return(se)
}

#' Inspect digest matches of a tximix-imported SE object
#'
#' @param se the SummarizedExperiment, or alternatively just
#' `metadata(se)$quantInfo`, a list of metadata
#' information from the quantification tool 
#' @param type what quantifier was used (see \code{\link{tximport}})
#' @param expanded_digest whether to include the full digest in the output, 
#' or just a shortened 6-char version
#' 
#' @return a tibble of the annotated and novel transcriptome information,
#' e.g. the index sequence digest, and if there is a match in the hash tables
#' 
#' @export
tximixInspectDigests <- function(se, type="oarfish", expanded_digest=FALSE) {
  
  # take from first sample
  if (is(se, "SummarizedExperiment")) {
    digestList <- metadata(se)$quantInfo$digest[,1]
  } else {
    # assume `se` isn't SE but the `quantInfo` item
    digestList <- se$digest[,1]
  }
  
  digests <- c(
    annotated = digestList$annotated_transcripts_digest$sha256_digests$sha256_seqs,
    novel = digestList$novel_transcripts_digest$sha256_digests$sha256_seqs
  )

  small_digest <- substr(digests, 1, 6)
 
  txomeInfo <- sapply(digests, getTxomeInfo, quiet=TRUE)

  out <- tibble(
    index=c("annotated","novel"), 
    source=NA, organism=NA, release=NA, 
    linkedTxome=NA, small_digest
  )
  if (expanded_digest) {
    out$digest <- digests
  }
  for (i in c("annotated","novel")) {
    if (!is.null(txomeInfo[[i]])) {
      cols <- c("source","organism","release","linkedTxome")
      out[match(i,out$index),cols] <- txomeInfo[[i]][cols]
    }
  }
  
  out
}

#' Update transcript metadatda for a tximix-imported SE object
#'
#' Will update the metadata on the SE object, using either
#' linkedTxome or linkedTxpData (preference to the former)
#'
#' @param se the SummarizedExperiment
#' @param txpData either GRanges or data.frame-type object
#' to use if there is not a match based on digest. 
#' This is used on a one-time basis, and transcripts
#' will be marked in metadata columns as `index = "user"``.
#' See `linkedTxome` or `linkedTxpData` for persistent
#' metadata storage/retrieval
#' @param ranges logical, whether to add rowRanges or rowData
#' @param order order in which to update the metadata
#' @param key the name of the column to use as the key
#' for merging metadata with the SE, which uses `rownames(se)`.
#' defaults to `tx_name` which often matches the transcript 
#' names in GENCODE
#'
#' @return a SummarizedExperiment with additional rowData,
#' or a RangedSummarizedExperiment
#'
#' @export
tximixUpdateTxpData <- function(
  se,
  txpData = NULL,
  ranges = FALSE,
  order = c("annotated", "novel"),
  key = "tx_name"
) {
  stopifnot(sort(order) == c("annotated", "novel"))

  # pull out digest list information from quantification tool
  digestList <- metadata(se)$quantInfo$digest[, 1]
  digests <- c(
    annotated = digestList$annotated_transcripts_digest$sha256_digests$sha256_seqs,
    novel = digestList$novel_transcripts_digest$sha256_digests$sha256_seqs
  )

  # pull out the txomeInfo for each
  txomeInfo <- sapply(digests, getTxomeInfo, quiet = TRUE)

  ranges_to_add <- GenomicRanges::GRanges()

  # in a specified order (default annotated then novel),
  # update the metadata columns in the rowData, which means
  # pulling out rowData, seeing what columns could be added/updated,
  # and then resaving to the rowData slot. This happens for each index.
  for (i in order) {
    if (!is.null(txomeInfo[[i]])) {
      # get the txdb
      txdb <- getTxDb(txomeInfo[[i]], useHub = FALSE, skipFtp = FALSE)
      txps <- getRanges(txdb = txdb, txomeInfo = txomeInfo[[i]], type = "txp")
      matches <- intersect(rownames(se), names(txps))
      if (length(matches) > 0) {
        message(paste0(
          "--",
          i,
          " index: adding metadata for ",
          length(matches),
          " transcripts"
        ))

        idx_txps <- match(matches, names(txps)) # index of the matches in the ranges
        txpDataToAdd <- mcols(txps)[idx_txps, ]

        # later in the function, ranges will be added
        if (ranges) {
          ranges_to_add <- c(ranges_to_add, txps[idx_txps])
        }

        # pull out rowData for metadata additions
        rowdata <- rowData(se)
        # if rowdata is totally empty, need to add one column
        if (ncol(rowdata) == 0) {
          rowdata[[key]] <- rownames(se)
        }
        # add in the metadata to the matching rows, and the index name
        rowdata <- mergeTxpDataIntoRowData(rowdata, txpDataToAdd, matches, indexName=i)

      } else {
        # matches of the transcripts from TxDb to the rows of SE was 0
        message(paste0("--", i, " index: no matching transcripts for the", ))
      }
    } else {
      # there was no linkedTxome to find
      message(
        paste0("--", i, " index: no transcript metadata found\n"),
        "  consider to add a 'linkedTxome', 'linkedTxpData'"
      )
    }
    # add the newly updated rowdata back to the SE
    SummarizedExperiment::rowData(se) <- rowdata
  }

  ### txpData ###

  # the above code looks up digests in the information stored with BiocFileCache, 
  # here user can provide 'txpData' on a one-time basis, labelled `index = "user"`
  if (!is.null(txpData)) {
    if (is(txpData, "GRanges")) {
      txps <- txpData # these will be used for the ranges
      names_txps <- names(txps) # used for matching
      txpDataToAdd <- mcols(txpData) # the metadata columns
      mcols(txps) <- NULL
    } else {
      # no ranges just data.frame like thing
      txpDataToAdd <- as(txpData, "DataFrame")
      names_txps <- txpData[[key]]
    }
    matches <- intersect(rownames(se), names_txps)
    if (length(matches) > 0) {
      message("txpData: adding transcript metadata for ", length(matches), " transcripts")
      idx_txps <- match(matches, names_txps) # index of the matches in the ranges
      txpDataToAdd <- txpDataToAdd[idx_txps, ] # put in order of matches
      if (ranges & is(txpData, "GRanges")) {
        ranges_to_add <- c(ranges_to_add, txps[idx_txps])
      }
      rowdata <- rowData(se)
      rowdata <- mergeTxpDataIntoRowData(rowdata, txpDataToAdd, matches, indexName="user")
      SummarizedExperiment::rowData(se) <- rowdata
    } else {
      if (is(txpData, "GRanges")) {
        message("txpData had no matching transcripts, check names(txpData)")
      } else {
        message("txpData had no matching transcripts, check 'key' column of txpData")
      }
      
    }
  }

  if (ranges) {
    # we've already dealt with metadata columns above, just add bare ranges
    mcols(ranges_to_add) <- NULL
    matches <- intersect(rownames(se), names(ranges_to_add))
    if (length(matches) < nrow(se)) {
      message(paste(
        "building RangedSE: subsetting to",
        length(matches),
        "out of",
        nrow(se),
        "rows with range data"
      ))
      ranges_to_add <- ranges_to_add[matches]
      se <- se[matches, ]
    }
    mcols(ranges_to_add) <- rowData(se)
    rowRanges(se) <- ranges_to_add
  }

  se
}

# txpDataToAdd and matches are in same order, not true for rowdata
mergeTxpDataIntoRowData <- function(rowdata, txpDataToAdd, matches, indexName) {
  # store the new transcript data back in the appropriate rows of the SE
  idx_rowdata <- match(matches, rownames(rowdata)) # index of the matches in the SE
  for (col in colnames(txpDataToAdd)) {
    if (!col %in% colnames(rowdata)) {
      rowdata[col] <- NA
    }
    rowdata[idx_rowdata, col] <- txpDataToAdd[, col]
  }
  # add the 'index' column and the indexName to the matching rows
  if (!"index" %in% colnames(rowdata)) {
    rowdata["index"] <- NA
  }
  rowdata[idx_rowdata, "index"] <- indexName
  rowdata
}
