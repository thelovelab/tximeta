#' tximeta
#'
#' tximeta is magic
#'
#' @param coldata a data.frame with at least two columns:
#' \itemize{
#' \item{"files"}{character vector pointing to sample quantification files}
#' \item{"names"}{character vector of sample names}
#' }
#' 
#' @return a SummarizedExperiment
#'
#' @importFrom tximport tximport
#' @importFrom rjson fromJSON
#' @importFrom AnnotationDbi loadDb saveDb
#'
#' @export
tximeta <- function(coldata) {

  stopifnot(all(c("files","names") %in% names(coldata)))

  files <- coldata$files
  names(files) <- coldata$names

  # get metadata from JSON files within quant dirs
  metaInfo <- lapply(files, getMetaInfo)
  indexSeqHash <- metaInfo[[1]]$index_seq_hash # first sample

  if (length(files) > 1) {
    hashes <- sapply(metaInfo, function(x) x$index_seq_hash)
    if (!all(hashes == indexSeqHash)) {
      stop("the samples do not share the same index, and cannot be imported")
    }
  }

  # try to import files early, so we don't waste user time
  # with metadata magic before a tximport error
  message("importing quantifications")
  txi <- tximport(files, type="salmon", txOut=TRUE)

  # now start building out metadata based on indexSeqHash
  hashtable <- read.csv(here("extdata","hashtable.csv"),stringsAsFactors=FALSE)
  idx <- match(indexSeqHash, hashtable$index_seq_hash)

  # did we find a match?
  # TODO: just return the SE instead of giving error here
  if (length(idx) == 0) stop("couldn't find metadata for the quant files")
  if (length(idx) > 1) stop("found more than one matching transcriptome...problem hash database")

  # now we can go get the GTF to annotate the ranges
  message(with(hashtable[idx,,drop=FALSE],
               paste0("found matching transcriptome:\n[ ",
                      source," - ",organism," - version ",version," ]")))  
  gtf <- hashtable$gtf[idx]
  txdbName <- paste0(basename(gtf),".sqlite")
  if (!file.exists(txdbName)) {
    message("building TxDb")
    txdb <- makeTxDbFromGFF(gtf)
    saveDb(txdb, file=txdbName)
  } else {
    message("loading existing TxDb")
    txdb <- loadDb(txdbName)
  }
  
  message("generating transcript ranges")
  txps <- transcripts(txdb)
  names(txps) <- txps$tx_name
  stopifnot(all(rownames(txi$abundance) %in% names(txps)))
  txps <- txps[rownames(txi$abundance)]
  
  message("fetching genome info")
  genome <- hashtable$genome[idx]
  chromInfo <- fetchExtendedChromInfoFromUCSC(genome2UCSC(genome))
  slstyle <- seqlevelsStyle(seqlevels(txps))
  slcol <- paste0(slstyle, "_seqlevel")
  chromInfoSub <- chromInfo[match(seqlevels(txps), chromInfo[[slcol]]),]
  seqlengths(seqinfo(txps)) <- chromInfoSub$UCSC_seqlength
  isCircular(seqinfo(txps)) <- chromInfoSub$circular
  genome(seqinfo(txps)) <- genome2UCSC(genome)

  # remove the files column from colData
  coldataSub <- subset(coldata, select=-files)

  # prepare metadata list
  tximetaInfo <- list(version=packageVersion("tximeta"),
                      importTime=Sys.time())

  txdbInfo <- metadata(txdb)$value
  names(txdbInfo) <- metadata(txdb)$name
  
  metaInfo <- reshapeMetaInfo(metaInfo)
  metadata <- list(
    quantInfo=metaInfo,
    txomeInfo=as.list(hashtable[idx,]),
    tximetaInfo=tximetaInfo,
    txdbInfo=txdbInfo
  )
  
  se <- SummarizedExperiment(assays=txi[c("abundance","counts","length")],
                             rowRanges=txps,
                             colData=coldataSub,
                             metadata=metadata)
  se
  
}

getMetaInfo <- function(file) {
  dir <- dirname(file)
  jsonPath <- file.path(dir, "aux_info/meta_info.json")
  stopifnot(file.exists(jsonPath))
  fromJSON(file=jsonPath)
}

reshapeMetaInfo <- function(metaInfo) {
  unionTags <- unique(unlist(lapply(metaInfo, names)))
  out <- lapply(unionTags, function(t) {
    sapply(seq_along(metaInfo), function(i) {
      metaInfo[[i]][[t]]
    })
  })
  names(out) <- unionTags
  if (all(out$eq_class_properties == list())) {
    out$eq_class_properties <- NULL
  }
  stopifnot(all(out$index_seq_hash == out$index_seq_hash[1]))
  stopifnot(all(out$index_name_hash == out$index_name_hash[1]))
  out$index_seq_hash <- out$index_seq_hash[1]
  out$index_name_hash <- out$index_name_hash[1]
  out
}

genome2UCSC <- function(x) {
  if (x == "GRCh38") {
    "hg38"
  } else {
    x
  }
}

#' Slightly less ugly liftOver functionality
#'
#' This includes the unlist and genome assignment
#' that needs to happen after liftOver
#'
#' @param ranges the incoming GRanges
#' @param chainfile a character vector pointing to a liftover chain file
#' @param to the name of the new genome
#'
#' @return a lifted GRanges
#'
#' @importFrom rtracklayer liftOver
#' 
#' @export
liftOverHelper <- function(ranges, chainfile, to) {
  chain <- import.chain(chainfile)
  out <- unlist(liftOver(ranges, chain))
  genome(out) <- to
  out
}
