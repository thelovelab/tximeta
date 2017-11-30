#' tximeta
#'
#' tximeta is magic
#'
#' @param coldata a data.frame with at least two columns:
#' \itemize{
#' \item{"files"}{character vector pointing to sample quantification files}
#' \item{"names"}{character vector of sample names}
#' }
#' @param ... initial arguments passed to \code{tximport}
#' 
#' @return a SummarizedExperiment
#'
#' @importFrom tximport tximport
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom AnnotationDbi loadDb saveDb
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcnew
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom rtracklayer import.chain liftOver
#'
#' @export
tximeta <- function(coldata, ...) {

  stopifnot(all(c("files","names") %in% names(coldata)))

  files <- coldata$files
  names(files) <- coldata$names

  # remove the files column from colData
  coldataSub <- subset(coldata, select=-files)

  # tximeta metadata
  tximetaInfo <- list(version=packageVersion("tximeta"),
                      importTime=Sys.time())
  
  # get quantifier metadata from JSON files within quant dirs
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
  txi <- tximport(files, type="salmon", txOut=TRUE, ...)

  # now start building out metadata based on indexSeqHash
  # TODO this is very temporary code obviously
  hashtable <- read.csv(here("extdata","hashtable.csv"),stringsAsFactors=FALSE)
  idx <- match(indexSeqHash, hashtable$index_seq_hash)

  metaInfo <- reshapeMetaInfo(metaInfo)
  metadata <- list(
    quantInfo=metaInfo,
    tximetaInfo=tximetaInfo
  )
  
  # did we find a match?
  if (is.na(idx)) {
    message("couldn't find matching transcriptome, returning un-ranged SummarizedExperiment") 
    se <- SummarizedExperiment(assays=txi[c("abundance","counts","length")],
                               colData=coldataSub,
                               metadata=metadata)
    return(se)
  }
  
  if (length(idx) > 1) stop("found more than one matching transcriptome...problem hash database")

  # now we can go get the GTF to annotate the ranges
  message(with(hashtable[idx,,drop=FALSE],
               paste0("found matching transcriptome:\n[ ",
                      source," - ",organism," - version ",version," ]")))  

  gtf <- hashtable$gtf[idx]
  txdbName <- basename(gtf)

  # TODO trial use of BiocFileCache to store the TxDb sqlite's
  bfc <- BiocFileCache(".")
  q <- bfcquery(bfc, txdbName)
  
  if (bfccount(q) == 0) {
    message("building TxDb")
    txdb <- makeTxDbFromGFF(gtf)
    savepath <- bfcnew(bfc, txdbName, ext="sqlite") 
    saveDb(txdb, file=savepath)
  } else {
    message(paste("loading existing TxDb created:",q$create_time[1]))
    txdb <- loadDb(bfcrpath(bfc, txdbName))
  }
  
  message("generating transcript ranges")
  txps <- transcripts(txdb)
  names(txps) <- txps$tx_name
  stopifnot(all(rownames(txi$abundance) %in% names(txps)))
  # TODO give a warning here if there are transcripts in TxDb not in Salmon index?
  txps <- txps[rownames(txi$abundance)]

  # TODO need a solution that doesn't rely on UCSC for the seqlevels
  # (this is not a good solution, just used for the prototype)
  # what is bad: the outgoing genome is now hg38 instead of GRCh38
  message("fetching genome info")
  genome <- hashtable$genome[idx]
  ucsc_genome <- genome2UCSC(genome)
  seqinfo(txps) <- Seqinfo(genome=ucsc_genome)[seqlevels(txps)]

  # add more metadata
  txdbInfo <- metadata(txdb)$value
  names(txdbInfo) <- metadata(txdb)$name
  metadata$txomeInfo <- as.list(hashtable[idx,])
  metadata$txdbInfo <- txdbInfo
  
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
  fromJSON(jsonPath)
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

# TODO obviously this will go
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
