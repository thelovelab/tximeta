#' Import quantification across mixed reference transcripts
#'
#' The _oarfish_ quantification tools allows a mix of 
#' `--annotated` reference transcripts (e.g. GENCODE, Ensembl) and 
#' `--novel` or custom transcripts (e.g. de novo assembled transcripts not present 
#' in the annotated set) to be used as the index for quantification.
#' `importData()` and associated functions facilitate import, reference identification, 
#' and addition of metadata across `annotated` and/or `novel` transcripts.
#' The `importData()` function alone imports the data, while inspection of the 
#' recognized digests and updating of transcript metadata is handled by subsequent functions
#' (listed in _See also_ section below).
#' 
#' oarfish with mixed reference transcript sets may have been generated with e.g.
#' 
#' \preformatted{
#' oarfish --only-index --annotated gencode.v48.transcripts.fa.gz \
#'   --novel my_novel_txps.fa.gz --seq-tech ont-cdna --threads 32 \
#'   --index-out gencode_plus_novel
#' oarfish --reads reads/experiment_rep1.fastq.gz --index gencode_plus_novel \
#'   --output quants/experiment_rep1 --seq-tech ont-cdna \
#'   --filter-group no-filters --threads 32
#' }
#'
#' @param coldata data.frame with columns `files` and `names` as in `tximeta()`
#' @param type what quantifier was used (see [tximport::tximport()]), for now 
#' `importData()` works for `"oarfish"` files
#' @param quiet whether to suppress printed messages
#' @param ... arguments passed to [tximport::tximport()]
#'
#' @return an un-ranged SummarizedExperiment (SE) object, for 
#' use with subsequent functions described in _See also_ section
#'
#' @seealso `inspectDigests()` and `updateMetadata()` for subsequent tasks
#' of inspecting digest matches and updating metadata, respectively.
#' `makeLinkedTxome()` can be used to add custom metadata into the registry used
#' for inspecting digests and then updating transcript data. A user may 
#' follow the workflow `importData()` > `inspectDigests()` > 
#' `makeLinkedTxome()` > `inspectDigests()` > `updateMetadata()`.
#' See also `makeLinkedTxpData()` for a lightweight alternative of linking
#' _GRanges_ metadata to a digest.
#' 
#' @examples
#' 
#' # oarfish files using a mix of --annotated and --novel transcripts
#' dir <- system.file("extdata/oarfish", package="tximportData")
#' names <- paste0("rep", 2:4)
#' files <- file.path(dir, paste0("sgnex_h9_", names, ".quant.gz"))
#' coldata <- data.frame(files, names)
#' 
#' # returns an un-ranged SE object
#' se <- importData(coldata, type="oarfish")
#' 
#' @export
importData <- function(coldata, type="oarfish", quiet=FALSE, ...) {

  if (!type == "oarfish") {
    warning("importData() supports oarfish files; broader support is planned in future updates")
  }
  
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

  digestList <- names(metaInfo[[1]]$digest)
  if (!all(paste0(c("annotated","novel"),"_transcripts_digest") %in% digestList))
      stop(
      "importData() is designed for mixed `annotated` and `novel` transcript references\n",
      "otherwise use tximeta() which will prioritize the `annotated` transcript set\n",
      "or tximeta(..., skipMeta=TRUE) to import all transcripts"
    )
  
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
    message("returning un-ranged SummarizedExperiment, see functions:\n",
            "-- inspectDigests() to check matching digests\n",
            "-- makeLinkedTxome/makeLinkedTxpData() to link digests to metadata\n",
            "-- updateMetadata() to update metadata and optionally add ranges"
          )

  return(se)
}

#' Inspect digest matches from `importData()` imported data
#' 
#' This function takes as input a _SummarizedExperiment_ as output by `importData()`
#' and returns a tibble with information about the digest-match status of two
#' indices (`annotated` and `novel`), with respect to _tximeta_ metadata.
#' Inspection of index digests can be run iteratively, checking if
#' the digests used in the mixed reference transcript set 
#' have a match against 1) pre-computed digests representing 
#' standard annotated sets (e.g. GENCODE, Ensembl, etc.) 
#' or 2) digests added by the user to a local registry with 
#' `makeLinkedTxome()` (GTF file)
#' or `makeLinkedTxpData` (_GRanges_-based metadata). 
#' Optional columns may be added if specified by 
#' `fullDigest=TRUE` (include the full digest) and/or 
#' `count=TRUE` (add matching transcript ID counts per index).
#' Following inspection, one can run `updateMetadata()` to automatically update
#' the transcript metadata using the sources indicated by this function.
#'
#' @param se the _SummarizedExperiment_ output by `importData()`,
#'  or alternatively just
#' `metadata(se)$quantInfo`, a list of metadata
#' information from the quantification tool 
#' (assuming `annotated` and `novel` indices both used)
#' @param type what quantifier was used (see [tximport::tximport()])
#' @param prefer vector of length up to 3, giving the preferred order of 
#' _tximeta_'s transcript registries to when finding matches, with elements:
#' `txome`: linkedTxome, 
#' `txpdata`: linkedTxpData,
#' `precomputed`: the pre-computed digests in tximeta
#' @param fullDigest whether to include the full digest string in the output, 
#' in addition to the shortened 6-char version
#' @param count whether to count the number of matching transcripts ID to each index
#' (only possible for those indices that have matching metadata).
#' Counting requires loading transcript data, either from locally cached databases
#' or from GTF files.
#' 
#' @return a 2-row tibble of the `annotated` and `novel` index, 
#' their matching information if available
#' (source, organism, release), for matches, 
#' whether it is a `linkedTxome` or a `linkedTxpData`
#' (both `FALSE`` for pre-computed) 
#' and a small 6 character version of the digest itself.
#' 
#' @examples
#' 
#' example(importData)
#' # now we have an `se` created by importData()...
#' inspectDigests(se)
#' # can then update the registry via makeLinkedTxome() and re-run inspection
#' 
#' @export
inspectDigests <- function(
  se, 
  type="oarfish", 
  prefer=c("txome","txpdata","precomputed"),
  fullDigest=FALSE, 
  count=FALSE
) {

  stopifnot(all(prefer %in% c("txome","txpdata","precomputed")))
  
  # take from first sample
  if (is(se, "SummarizedExperiment")) {
    digestList <- metadata(se)$quantInfo$digest[,1]
  } else {
    stopifnot(!count) # counting transcripts to indices requires rownames of an SE
    # assume `se` isn't SE but the `quantInfo` item
    digestList <- se$digest[,1]
  }
  
  # need to check, even though importData would have thrown error
  stopifnot(all(paste0(c("annotated","novel"),"_transcripts_digest") %in% names(digestList)))

  digests <- c(
    annotated = digestList$annotated_transcripts_digest$sha256_digests$sha256_seqs,
    novel = digestList$novel_transcripts_digest$sha256_digests$sha256_seqs
  )

  small_digest <- substr(digests, 1, 6)
 
  txomeInfo <- lapply(digests, getTxomeInfo, prefer, quiet=TRUE)

  # this is the tibble the function will return
  out <- tibble(
    index=c("annotated","novel"), 
    source=NA, organism=NA, release=NA, genome=NA,
    linkedTxome=NA, linkedTxpData=NA, small_digest
  )

  # put in the full digest if requested
  if (fullDigest) {
    out$digest <- digests
  }

  # columns to pull from the txomeInfo item
  cols <- c("source","organism","release","genome","linkedTxome","linkedTxpData")
  for (i in c("annotated","novel")) {
    # if there is a txomeInfo match, populate the outgoing tibble
    if (!is.null(txomeInfo[[i]])) {
      out[match(i,out$index),cols] <- txomeInfo[[i]][cols]
    }
  }

  if (count) {
    out$count <- 0
    for (i in c("annotated","novel")) {
      if (!is.null(txomeInfo[[i]])) {
        txps <- getTxpsFromTxome(txomeInfo = txomeInfo[[i]])
        out[match(i,out$index),"count"] <- sum(names(txps) %in% rownames(se))
      }
    }
  }

  out
}

#' Update transcript metadatda for `importData()` imported data
#'
#' This function takes as input a _SummarizedExperiment_ as output by `importData()`,
#' and will update the metadata on the transcripts when possible 
#' (updating `rowData` and/or `rowRanges` depending on the value of `ranges`).
#' `importData()` uses metadata pulled from digest matches in registries used by _tximeta_
#' (`linkedTxome`, `linkedTxpData`, and the pre-computed digests).
#' Additionally, _GRanges_ or _data.frame_-type data can be provided on a one-time basis 
#' via the argument `txpData`, which will annotate transcripts with `index="user"`.
#' See `inspectDigests()` for how to inspect which indices have matching digests, 
#' and how to link data to local metadata in a persistent manner.
#'
#' @param se the _SummarizedExperiment_ (SE) output by `importData()`
#' @param txpData either _GRanges_ or _data.frame_-type object
#' to use if there is not a match based on digest. 
#' This is used on a one-time basis, and transcripts
#' will be marked in metadata columns as `index = "user"``.
#' See `makeLinkedTxome()` or `makeLinkedTxpData()` for persistent
#' metadata storage/retrieval
#' @param ranges logical, whether to add `rowRanges` (or just `rowData`)
#' @param prefer vector of length up to 3, giving the preferred order of 
#' _tximeta_'s transcript registries to when finding matches, with elements:
#' `txome`: linkedTxome, 
#' `txpdata`: linkedTxpData,
#' `precomputed`: the pre-computed digests in tximeta
#' @param order order of index, in which to update the metadata, by default 
#' the order is `annotation`, then `novel`, then `user`, info supplied 
#' here as `txpData`
#' @param key a named character vector of length 3. For each index
#' (annotated, novel, and user) `key` is the name of the column to 
#' use for merging metadata with `rownames(se)`.
#' The `user` index corresponds to data provided here as `txpData`
#' Defaults to `"tx_name"` which often matches the transcript 
#' names in GENCODE
#'
#' @return a _SummarizedExperiment_ with new `rowData`,
#' or a _RangedSummarizedExperiment_ with new metadata
#'
#' @examples
#' 
#' example(importData)
#' 
#' # build custom novel GRanges data
#' library(GenomicRanges)
#' novel <- data.frame(
#'   seqnames = paste0("chr", rep(1:22, each=500)),  
#'   start = 1e6 + 1 + 0:499 * 1000, end = 1e6 + 1 + 0:499 * 1000 + 1000 - 1,
#'   strand = "+", tx_name = paste0("novel", 1:(22*500)),
#'   gene_id = paste0("novel_gene", rep(1:(22*10), each=50)), type = "protein_coding"
#' )
#' novel_gr <- as(novel, "GRanges")
#' names(novel_gr) <- novel$tx_name
#' 
#' # now update the metadata + ranges:
#' \dontrun{
#' # this requires connection to internet (will download GENCODE GTF via FTP)
#' se_with_ranges <- updateMetadata(
#'   se, txpData=novel_gr, ranges=TRUE
#' )
#' mcols(se_with_ranges)
#' }
#' 
#' @export
updateMetadata <- function(
  se,
  txpData = NULL,
  ranges = FALSE,
  prefer=c("txome","txpdata","precomputed"),
  order = c("annotated", "novel", "user"),
  key = c(annotated="tx_name", novel="tx_name", user="tx_name")
) {

  # check our standard index names
  idx_nms <- c("annotated","novel")
  all_idx_nms <- c(idx_nms, "user")
  stopifnot(all(order %in% all_idx_nms))
  stopifnot(all(names(key) %in% all_idx_nms))

  # pull out our digest's list of info:
  # this lives in metadata as information coming from the quantification tool
  digestList <- metadata(se)$quantInfo$digest[, 1]
  stopifnot(all(paste0(c("annotated","novel"),"_transcripts_digest") %in% names(digestList)))

  # we will just use a named list for `digest` here 
  # (maybe taking advantage of partial matching from seqcol later)
  digests <- c(
    annotated = digestList$annotated_transcripts_digest$sha256_digests$sha256_seqs,
    novel = digestList$novel_transcripts_digest$sha256_digests$sha256_seqs
  )

  # pull out the txomeInfo for each index
  txomeInfo <- lapply(digests, getTxomeInfo, prefer=c("txome","txpdata","precomputed"), quiet = TRUE)

  # empty GRanges, add to this per index / txpData in loop below
  if (ranges) {
    ranges_to_add <- GenomicRanges::GRanges()
  }

  # Proceed for each index `i` in a user-specified `order`, 
  # (by default annotated then novel then `txpData`),
  # updating the metadata columns in the rowData:
  #  
  # - pulling out rowData, 
  # - seeing what columns could be added/updated,
  # - resaving to the rowData slot. 
  # 
  # then move to the next index.
  for (i in order) {

    # get transcript data and determine 
    # the matches of rownames of `se` to these
    matches <- c()
    # first, annotated / novel index
    if (i %in% idx_nms) {
      if (!is.null(txomeInfo[[i]])) {
        txps <- getTxpsFromTxome(txomeInfo = txomeInfo[[i]])
        names_txps <- names(txps) # used for matching later
        txpDataToAdd <- mcols(txps) # metadata columns to work with
        matches <- intersect(rownames(se), names_txps)
      } else {
        # no match
        message(
          paste0("--", i, ": no transcript metadata found\n"),
          "  consider using `linkedTxome`, or `linkedTxpData` (see man pages)"
        )
      }

      # for the txpData-provided information..
    } else if (i == "user") {
      # either `txpData` is GRanges or data.frame-like thing
      if (!is.null(txpData)) {
        if (is(txpData, "GRanges")) {
          txps <- txpData # these will be used for the ranges
          names_txps <- names(txps) # used for matching later
          txpDataToAdd <- mcols(txpData) # metadata columns to work with
          mcols(txps) <- NULL
        } else {
          # no ranges just data.frame-like thing
          txpDataToAdd <- as(txpData, "DataFrame")
          names_txps <- txpData[[key[i]]]
        }
        matches <- intersect(rownames(se), names_txps)
      } else {

        # txpData wasn't provided, we are in the "user" part of the loop
        # need to skip to the end... and not do the next part again
        break

      }
    }

    # now we have `txps`/`txpDataToAdd` and `matches`
    if (length(matches) > 0) {
      message(paste0(
        "--",
        i,
        ": adding metadata for ",
        length(matches),
        " transcripts"
      ))
      idx_txps <- match(matches, names_txps) # index of the matches in the ranges
      txpDataToAdd <- txpDataToAdd[idx_txps, ] # put in order of matches
      
      # later in the function, ranges will be added
      if (ranges) {
        if (i %in% idx_nms | is(txpData, "GRanges")) {
          ranges_to_add <- c(ranges_to_add, txps[idx_txps])
        }
      }

      # pull out rowData for metadata additions
      rowdata <- rowData(se)
      # if rowdata is totally empty, need to add one column
      if (ncol(rowdata) == 0) {
        rowdata[[key[i]]] <- rownames(se)
      }
      # add in the metadata to the matching rows, and the index name
      rowdata <- mergeTxpDataIntoRowData(rowdata, txpDataToAdd, matches, indexName=i)
    } else {
      # matches of the transcripts from various sources to the rows of SE was 0
      message(paste0("--", i, ": no matching transcripts"))
    }

    # finally add the newly updated rowdata back to the SE
    SummarizedExperiment::rowData(se) <- rowdata

    # go to next index...
  }

  if (ranges) {
    # we've already dealt with metadata columns above, just add bare ranges
    mcols(ranges_to_add) <- NULL
    rng_matches <- intersect(rownames(se), names(ranges_to_add))
    if (length(rng_matches) < nrow(se)) {
      # print a note that we are subsetting to a smaller RangedSE than input
      message(paste(
        "building RangedSE: subsetting to",
        length(rng_matches),
        "out of",
        nrow(se),
        "rows with range data"
      ))
    }
    ranges_to_add <- ranges_to_add[rng_matches]
    se <- se[rng_matches, ]
    mcols(ranges_to_add) <- rowData(se)
    rowRanges(se) <- ranges_to_add
  }

  se
}

### un-exported helper functions ###

# txpDataToAdd and matches are in same order, not true for rowdata
mergeTxpDataIntoRowData <- function(rowdata, txpDataToAdd, matches, indexName) {
  # store the new transcript data back in the appropriate rows of the SE
  idx_rowdata <- match(matches, rownames(rowdata)) # index of the matches in the SE
  for (col in colnames(txpDataToAdd)) {
    if (!col %in% colnames(rowdata)) {
      # initialize with NA
      vector <- txpDataToAdd[, col]
      vector <- S4Vectors::endoapply(vector, \(x) NA)
      rowdata[col] <- rep(vector, length.out = nrow(rowdata))
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

# function to pull out GRanges txps from txomeInfo whether 
# linkedTxome, pre-computed, or linkedTxpData
getTxpsFromTxome <- function(txomeInfo) {
  if (!txomeInfo$linkedTxpData) {
    # we have a linkedTxome or pre-computed digest match
    suppressMessages({
      txdb <- getTxDb(txomeInfo, useHub = FALSE, skipFtp = FALSE)
      txps <- getRanges(txdb = txdb, txomeInfo = txomeInfo, type = "txp")
    })
  } else if (txomeInfo$linkedTxpData) {
    # we have a linkedTxpData match
    digest32 <- substr(txomeInfo$digest, 1, 32)
    txpDataName <- paste0("txpdata-", digest32)
    bfc <- BiocFileCache(getBFCLoc())
    if (!existsInBFC(txpDataName, bfc)) {
      stop(paste0("TxpData of name: [", txpDataName, "] was expected in BFC"))
    }
    loadpath <- bfcrpath(bfc, rnames = txpDataName)
    txps <- readRDS(loadpath)
  }
  txps
}
