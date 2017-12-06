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
#' @importFrom tximport tximport summarizeToGene
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom AnnotationDbi loadDb saveDb
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts
#' @importFrom ensembldb ensDbFromGtf EnsDb
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcnew bfccount bfcrpath
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
  # reshape
  metaInfo <- reshapeMetaInfo(metaInfo)
  # start to build metadata list
  metadata <- list(
    quantInfo=metaInfo,
    tximetaInfo=tximetaInfo
  )
  # try to import files early, so we don't waste user time
  # with metadata magic before a tximport error
  message("importing quantifications")
  txi <- tximport(files, type="salmon", txOut=TRUE, ...)

  # TODO need to store "countsFromAbundance" value in outgoing 'se'

  # try and find a matching txome
  txomeInfo <- getTxomeInfo(indexSeqHash)
  if (is.null(txomeInfo)) {
    message("couldn't find matching transcriptome, returning un-ranged SummarizedExperiment") 
    se <- SummarizedExperiment(assays=txi[c("counts","abundance","length")],
                               colData=coldataSub,
                               metadata=metadata)
    return(se)
  }

  # go build or load a TxDb from the gtf
  txdb <- getTxDb(txomeInfo)
  
  message("generating transcript ranges")
  txps <- transcripts(txdb)
  names(txps) <- txps$tx_name
  stopifnot(all(rownames(txi$abundance) %in% names(txps)))
  
  # TODO give a warning here if there are transcripts in TxDb not in Salmon index?
  # ...hmm, maybe not because that is now a "feature" given derivedTxomes that
  # are a simple subset of the transcripts/genes in the source FASTA & GTF
  txps <- txps[rownames(txi$abundance)]

  # Ensembl already has nice seqinfo attached, if not:
  if (txomeInfo$source != "Ensembl") {
    # TODO can we get a solution that doesn't rely on UCSC for the seqlevels?
    message("fetching genome info")
    ucsc_genome <- genome2UCSC(txomeInfo$genome)
    seqinfo(txps) <- Seqinfo(genome=ucsc_genome)[seqlevels(txps)]
  }

  # add more metadata
  txdbInfo <- metadata(txdb)$value
  names(txdbInfo) <- metadata(txdb)$name
  metadata$txomeInfo <- txomeInfo
  metadata$txdbInfo <- txdbInfo
  
  se <- SummarizedExperiment(assays=txi[c("counts","abundance","length")],
                             rowRanges=txps,
                             colData=coldataSub,
                             metadata=metadata)
  se
  
}

getMetaInfo <- function(file) {
  dir <- dirname(file)
  jsonPath <- file.path(dir,"aux_info","meta_info.json")
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

getTxomeInfo <- function(indexSeqHash) {
  # now start building out metadata based on indexSeqHash
  # TODO this is very temporary code obviously
  hashtable <- read.csv(here("extdata","hashtable.csv"),stringsAsFactors=FALSE)
  m <- match(indexSeqHash, hashtable$index_seq_hash)
  if (!is.na(m)) {
    # now we can go get the GTF to annotate the ranges
    txomeInfo <- as.list(hashtable[m,])
    message(with(txomeInfo,
                 paste0("found matching transcriptome:\n[ ",
                        source," - ",organism," - version ",version," ]")))
    return(txomeInfo)
  } 
  # TODO trial use of BiocFileCache to store the TxDb sqlite's
  # later change this to the default BiocFileCache(), or figure
  # out how users can point to their own preferred BFC locations
  bfc <- BiocFileCache(".")
  q <- bfcquery(bfc, "derivedTxomeDF")
  stopifnot(bfccount(q) < 2)
  if (bfccount(q) == 1) {
    loadpath <- bfcrpath(bfc, "derivedTxomeDF")
    derivedTxomeDF <- readRDS(loadpath)
    m2 <- match(indexSeqHash, derivedTxomeDF$index_seq_hash)
    if (!is.na(m2)) {
      txomeInfo <- as.list(derivedTxomeDF[m2,])
      message(with(txomeInfo,
                   paste0("found matching derived transcriptome:\n[ ",
                          source," - ",organism," - version ",version," ]")))    
      return(txomeInfo)
    }
  }
  return(NULL)
}

getTxDb <- function(txomeInfo) {
  txdbName <- basename(txomeInfo$gtf)  
  bfc <- BiocFileCache(".")
  q <- bfcquery(bfc, txdbName)  
  if (bfccount(q) == 0) {
    savepath <- bfcnew(bfc, txdbName, ext="sqlite") 
    if (txomeInfo$source == "Ensembl") {
      message("building EnsDb with 'ensembldb' package")
      # TODO suppress warnings here or what?
      suppressWarnings(ensDbFromGtf(txomeInfo$gtf, outfile=savepath))
      txdb <- EnsDb(savepath)
    } else {
      message("building TxDb with 'GenomicFeatures' package")
      txdb <- makeTxDbFromGFF(txomeInfo$gtf)
      saveDb(txdb, file=savepath)
    }
  } else {
    loadpath <- bfcrpath(bfc, txdbName)
    if (txomeInfo$source == "Ensembl") {
      message(paste("loading existing EnsDb created:",q$create_time[1]))
      txdb <- EnsDb(loadpath)
    } else {
      message(paste("loading existing TxDb created:",q$create_time[1]))
      txdb <- loadDb(loadpath)
    }
  }
  txdb
}

