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
#' @param ... passed to tximport
#'
#' @return an unranged SummarizedExperiment
#'
#' @export
tximix <- function(coldata, type="oarfish", ...) {
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
  
  return(se)
}

#' Inspect the annotations of a tximix-imported SE object
#'
#' @param se the SummarizedExperiment
#' @param type what quantifier was used (see \code{\link{tximport}})
#'
#' @return a tibble of the annotated and novel transcriptome information,
#' e.g. the hash, and if there is a match in the hash tables
#' 
#' @export
tximixInspect <- function(se, type="oarfish") {
  
  # take from first sample
  indexSeqHashList <- metadata(se)$quantInfo$digest[,1]

  indexSeqHashes <- c(
    indexSeqHashList$annotated_transcripts_digest$sha256_digests$sha256_seqs,
    indexSeqHashList$novel_transcripts_digest$sha256_digests$sha256_seqs
  )

  indexSeqSubstr <- substr(indexSeqHashes, 1, 6)

  tibble(index=c("annotated","novel"), indexSeqSubstr, indexSeqHashes)
  
}
