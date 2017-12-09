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
#' @importFrom rappdirs user_cache_dir
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
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
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
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
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

#' @rdname getTximetaBFC
#' 
#' @export
setTximetaBFC <- function() {
  message("Which BiocFileCache directory should tximeta use? (press Enter to cancel)")
  bfcloc <- file.choose()
  tximetaDir <- user_cache_dir("tximeta")
  bfclocFile <- file.path(tximetaDir, "bfcloc.txt")
  if (!file.exists(tximetaDir)) dir.create(tximetaDir)
  write(bfcloc, file=bfclocFile)
  bfcloc
}

#' Get or set the directory of the BiocFileCache used by tximeta
#'
#' Running \code{getTximetaBFC} will report the saved directory,
#' if it has been determined, or will return NULL.
#' Running \code{setTximetaBFC} will ask the user to specify a
#' BiocFileCache directory for accessing and saving TxDb sqlite files.
#'
#' @return the directory of the BiocFileCache used by tximeta
#'
#' @rdname getTximetaBFC
#' 
#' @export
getTximetaBFC <- function() {
  tximetaDir <- user_cache_dir("tximeta")
  bfclocFile <- file.path(tximetaDir, "bfcloc.txt")
  if (!file.exists(bfclocFile)) {
    message("tximeta's BiocFileCache location has not yet been set")
    NULL
  } else {
    scan(bfclocFile, what="char", quiet=TRUE)
  }
}

getBFCLoc <- function() {
  defaultDir <- user_cache_dir(appname="BiocFileCache")
  prompt <- paste("",
  "tximeta needs a BiocFileCache directory to access and save TxDb objects.",
  paste0("Do you wish to use the default directory: '",defaultDir,"'?"),
  "If not, a temporary directory that is specific to this R session will be used.","",
  "You can always change this directory later by running: setTximetaBFC()",
  "Or enter [0] to exit and set this directory manually now.",
  sep="\n")

  # this file tells us which BFC dir has been previously chosen use with tximeta
  tximetaDir <- user_cache_dir("tximeta")
  bfclocFile <- file.path(tximetaDir, "bfcloc.txt")
  if (!file.exists(bfclocFile)) {
    if (interactive()) {
      ans <- menu(c("Yes (use default)", "No (use temp)"), title=prompt)
      if (ans == 0) stop("no BiocFileCache directory choice made at this time")
      if (ans == 1) {
        bfcloc <- defaultDir
      } else {
        bfcloc <- tempdir()
      }
      if (!file.exists(tximetaDir)) dir.create(tximetaDir)
      write(bfcloc, file=bfclocFile)
    } else {
      bfcloc <- tempdir()
    }
  } else {
    bfcloc <- scan(bfclocFile, what="char", quiet=TRUE)
  }
  bfcloc
}
