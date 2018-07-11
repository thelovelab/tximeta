#' tximeta: Import transcript abundances with automagic population of metadata
#'
#' \code{tximeta} leverages the hash signature of the Salmon or Sailfish index,
#' in addition to a number of core Bioconductor packages (GenomicFeatures,
#' ensembldb, GenomeInfoDb, BiocFileCache) to automatically populate metadata
#' for the user, without additional effort from the user. Note that this package
#' is in "beta" / under development.
#'
#' \code{tximeta} checks the hash signature of the index against a database
#' of known transcriptomes (this database under construction) or a locally stored
#' \code{derivedTxome} (see \code{link{makeDerivedTxome}}), and then will
#' automatically populate, e.g. the transcript locations, the transcriptome version,
#' the genome with correct chromosome lengths, etc. It allows for automatic
#' and correct summarization of transcript-level quantifications to the gene-level
#' via \code{\link{summarizeToGene}} without the need to manually build
#' a \code{tx2gene} table.
#'
#' \code{tximeta} on the first run will ask where the BiocFileCache for
#' this package should be kept, either using a default location or a temporary
#' directory. At any point, the user can specify a location using
#' \code{\link{setTximetaBFC}} and this choice will be saved for future sessions.
#' Multiple users can point to the same BiocFileCache, such that
#' transcript databases (TxDb) associated with certain Salmon or Sailfish indices
#' and \code{derivedTxomes} can be accessed by different users without additional
#' effort or time spent downloading/building the relevant TxDb.
#' If the BiocFileCache location is not the default one in the user directory,
#' then \code{tximeta} will set the permissions to \code{666} so that
#' multiple users can read and write.
#'
#' @param coldata a data.frame with at least two columns:
#' \itemize{
#' \item{"files"}{character vector pointing to sample quantification files}
#' \item{"names"}{character vector of sample names}
#' }
#' @param ... initial arguments passed to \code{tximport}
#' 
#' @return a SummarizedExperiment with metadata on the \code{rowRanges}.
#' (if the hash signature in the Salmon or Sailfish index does not match
#' any known transcriptomes, or any locally saved \code{derivedTxome},
#' \code{tximeta} will just return a non-ranged SummarizedExperiment)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays colData
#' @importFrom S4Vectors metadata mcols mcols<-
#' @importFrom tximport tximport summarizeToGene
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom AnnotationDbi loadDb saveDb select keys mapIds
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts genes
#' @importFrom ensembldb ensDbFromGtf EnsDb
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcnew bfccount bfcrpath
#' @importFrom GenomeInfoDb Seqinfo genome<- seqinfo<- seqlevels
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom rappdirs user_cache_dir
#' @importFrom utils menu packageVersion read.csv
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
  # best this would be an external data package
  hashfile <- file.path(system.file("extdata",package="tximeta"),"hashtable.csv")
  hashtable <- read.csv(hashfile,stringsAsFactors=FALSE)
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
  # want to enable multiple users to read/write to the tximeta BiocFileCache
  # this should perhaps only happen when the bfloc is not the user-based location
  if (bfcloc != user_cache_dir(appname="BiocFileCache")) {
    bfcSqlite <- file.path(bfcloc, "BiocFileCache.sqlite")
    if (file.info(bfcSqlite)$mode != "666") {
      Sys.chmod(paths=bfcSqlite, mode="666", use_umask=FALSE)
    }
  }
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

