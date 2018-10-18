#' tximeta: Transcript quantification import with automatic metadata
#'
#' \code{tximeta} leverages the hash signature of the Salmon or Sailfish index,
#' in addition to a number of core Bioconductor packages (GenomicFeatures,
#' ensembldb, GenomeInfoDb, BiocFileCache) to automatically populate metadata
#' for the user, without additional effort from the user. Note that this package
#' is in "beta" / under development.
#'
#' Most of the code in \code{tximeta} works to add metadata and transcript ranges
#' when the quantification was performed with Salmon or Sailfish. However,
#' \code{tximeta} can be used with any quantification \code{type} that is supported
#' by \code{\link{tximport}}, where it will return an un-ranged SummarizedExperiment.
#' 
#' \code{tximeta} checks the hash signature of the index against a database
#' of known transcriptomes (this database under construction) or a locally stored
#' \code{linkedTxome} (see \code{link{makeLinkedTxome}}), and then will
#' automatically populate, e.g. the transcript locations, the transcriptome release,
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
#' and \code{linkedTxomes} can be accessed by different users without additional
#' effort or time spent downloading/building the relevant TxDb.
#'
#' In order to allow that multiple users can read and write to the
#' same location, one should set the BiocFileCache directory to
#' have group write permissions (g+w).
#'
#' @param coldata a data.frame with at least two columns (others will propogate to object):
#' \itemize{
#' \item{\code{files} - character, paths of quantification files}
#' \item{\code{names} - character, sample names}
#' }
#' if \code{coldata} is a vector, it is assumed to be the paths of quantification files
#' and unique sample names are created
#' @param type what quantifier was used (see \code{\link{tximport}})
#' @param ... arguments passed to \code{tximport}
#' 
#' @return a SummarizedExperiment with metadata on the \code{rowRanges}.
#' (if the hash signature in the Salmon or Sailfish index does not match
#' any known transcriptomes, or any locally saved \code{linkedTxome},
#' \code{tximeta} will just return a non-ranged SummarizedExperiment)
#'
#' @examples
#'
#' # point to a Salmon quantification file:
#' dir <- system.file("extdata/salmon_dm", package="tximportData")
#' files <- file.path(dir, "SRR1197474_cdna", "quant.sf.gz") 
#' coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
#'
#' # normally we would just run the following which would download the appropriate metadata
#' # se <- tximeta(coldata)
#'
#' # for this example, we instead point to a local path where the GTF can be found
#' # by making a linkedTxome:
#' dir <- system.file("extdata", package="tximeta")
#' indexDir <- file.path(dir, "Drosophila_melanogaster.BDGP6.cdna.v92_salmon_0.10.2")
#' fastaFTP <- "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
#' dir2 <- system.file("extdata/salmon_dm", package="tximportData")
#' gtfPath <- file.path(dir2,"Drosophila_melanogaster.BDGP6.92.gtf.gz")
#' makeLinkedTxome(indexDir=indexDir, source="Ensembl", organism="Drosophila melanogaster",
#'                 release="92", genome="BDGP6", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#' se <- tximeta(coldata)
#'
#' # to clear the entire linkedTxome table
#' # (don't run unless you want to clear this table!)
#' # bfcloc <- getTximetaBFC()
#' # bfc <- BiocFileCache(bfcloc)
#' # bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays assayNames colData
#' @importFrom S4Vectors metadata mcols mcols<-
#' @importFrom tximport tximport summarizeToGene
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom AnnotationDbi loadDb saveDb select keys mapIds
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts genes
#' @importFrom ensembldb ensDbFromGtf EnsDb
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcnew bfccount bfcrpath
#' @importFrom tibble tibble
#' @importFrom GenomeInfoDb Seqinfo genome<- seqinfo<- seqlevels
#' @importFrom rappdirs user_cache_dir
#' @importFrom utils menu packageVersion read.csv
#' @importFrom methods is
#'
#' @export
tximeta <- function(coldata, type="salmon", ...) {

  if (is(coldata, "vector")) {
    coldata <- data.frame(files=coldata, names=seq_along(coldata))
  }
  
  stopifnot(all(c("files","names") %in% names(coldata)))
  
  files <- as.character(coldata$files)
  names(files) <- coldata$names

  if (!all(file.exists(files))) {
    stop("the files do not exist at the location specified by 'coldata$files'")
  }
  
  # remove the files column from colData
  coldata <- subset(coldata, select=-files)

  # tximeta metadata
  tximetaInfo <- list(version=packageVersion("tximeta"),
                      importTime=Sys.time())

  metadata <- list(tximetaInfo=tximetaInfo)
  
  if (!type %in% c("salmon","sailfish")) {
    txi <- tximport(files, type=type, ...)
    metadata$countsFromAbundance <- txi$countsFromAbundance
    se <- makeUnrangedSE(txi, coldata, metadata)
    return(se)
  }
  
  # get quantifier metadata from JSON files within quant dirs
  metaInfo <- lapply(files, getMetaInfo)
  # Salmon's SHA-256 hash of the index is called "index_seq_hash" in the meta_info.json file
  indexSeqHash <- metaInfo[[1]]$index_seq_hash # first sample  
  if (length(files) > 1) {
    hashes <- sapply(metaInfo, function(x) x$index_seq_hash)
    if (!all(hashes == indexSeqHash)) {
      stop("the samples do not share the same index, and cannot be imported")
    }
  }
  # reshape
  metaInfo <- reshapeMetaInfo(metaInfo)
  # add to metadata list
  metadata$quantInfo <- metaInfo
  
  # try to import files early, so we don't waste user time
  # with metadata magic before a tximport error
  message("importing quantifications")
  txi <- tximport(files, type=type, txOut=TRUE, ...)
  metadata$countsFromAbundance <- txi$countsFromAbundance

  # try and find a matching txome
  txomeInfo <- getTxomeInfo(indexSeqHash)
  if (is.null(txomeInfo)) {
    message("couldn't find matching transcriptome, returning un-ranged SummarizedExperiment")
    se <- makeUnrangedSE(txi, coldata, metadata)
    return(se)
  }

  # go build or load a TxDb from the gtf
  txdb <- getTxDb(txomeInfo)
  
  message("generating transcript ranges")
  # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?
  if (txomeInfo$source == "Ensembl") {
    suppressWarnings({
        txps <- transcripts(txdb)
    })
  } else {
    suppressWarnings({
        txps <- transcripts(txdb, columns=c("tx_id","gene_id","tx_name"))
    })
  }
  names(txps) <- txps$tx_name

  # put 'counts' in front to facilitate DESeqDataSet construction
  assays <- txi[c("counts","abundance","length")]

  # if there are inferential replicates or inferential variance
  if ("infReps" %in% names(txi)) {
    infReps <- rearrangeInfReps(txi)
    infReps <- lapply(infReps, function(mat) {
      rownames(mat) <- rownames(assays[["counts"]])
      colnames(mat) <- colnames(assays[["counts"]])
      mat
    })
    assays <- c(assays, infReps)
  } else if ("variance" %in% names(txi)) {
    assays <- c(assays, txi["variance"])
  }
  
  # TODO temporary hack: Ensembl FASTA has txp version, Ensembl GTF it is not in the txname
  # meanwhile Gencode GTF puts the version in the name
  if (txomeInfo$source == "Ensembl") {
    txId <- sub("\\..*", "", rownames(assays[["counts"]]))
    for (nm in names(assays)) {
      rownames(assays[[nm]]) <- txId
    }
  }

  assays <- checkAssays2Txps(assays, txps)
  
  # TODO give a warning here if there are transcripts in TxDb not in Salmon index?
  # ...hmm, maybe not because that is now a "feature" given linkedTxomes that
  # are a simple subset of the transcripts/genes in the source FASTA & GTF
  txps <- txps[rownames(assays[["counts"]])]

  # Ensembl already has nice seqinfo attached, if not:
  if (txomeInfo$source != "Ensembl") {
    # TODO can we get a solution that doesn't rely on UCSC for the seqlevels?
    # this produces an error if not connected to internet
    message("fetching genome info")
    ucsc_genome <- genome2UCSC(txomeInfo$genome)
    seqinfo(txps) <- Seqinfo(genome=ucsc_genome)[seqlevels(txps)]
  }

  # add more metadata
  txdbInfo <- metadata(txdb)$value
  names(txdbInfo) <- metadata(txdb)$name
  metadata$txomeInfo <- txomeInfo
  metadata$txdbInfo <- txdbInfo

  se <- SummarizedExperiment(assays=assays,
                             rowRanges=txps,
                             colData=coldata,
                             metadata=metadata)
  se
  
}

# read metadata files from Salmon directory
getMetaInfo <- function(file) {
  dir <- dirname(file)
  jsonPath <- file.path(dir,"aux_info","meta_info.json")
  stopifnot(file.exists(jsonPath))
  fromJSON(jsonPath)
}

# reshape metadata info from Salmon
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

# temporary function to map from GRCh38 to hg38 to allow easy
# comparison with UCSC objects from AnnotationHub...
# TODO obviously this will have to go/be rethought
genome2UCSC <- function(x) {
  if (x == "GRCh38") {
    "hg38"
  } else {
    x
  }
}

# build out metadata based on indexSeqHash
getTxomeInfo <- function(indexSeqHash) {

  # first try to find any linkedTxomes in the linkedTxomeTbl
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, "linkedTxomeTbl")
  # there should only be one such entry in the tximeta bfc
  stopifnot(bfccount(q) < 2)
  if (bfccount(q) == 1) {
    loadpath <- bfcrpath(bfc, "linkedTxomeTbl")
    linkedTxomeTbl <- readRDS(loadpath)
    m <- match(indexSeqHash, linkedTxomeTbl$sha256)
    if (!is.na(m)) {
      txomeInfo <- as.list(linkedTxomeTbl[m,])
      message(paste0("found matching linked transcriptome:\n[ ",
                     txomeInfo$source," - ",txomeInfo$organism," - release ",txomeInfo$release," ]"))
      return(txomeInfo)
    }
  }

  # if not in linkedTxomes try the pre-computed hashtable...

  # TODO this is very temporary code obviously
  # best this would be an external data package
  hashfile <- file.path(system.file("extdata",package="tximeta"),"hashtable.csv")
  hashtable <- read.csv(hashfile,stringsAsFactors=FALSE)
  m <- match(indexSeqHash, hashtable$sha256)
  if (!is.na(m)) {
    # now we can go get the GTF to annotate the ranges
    txomeInfo <- as.list(hashtable[m,])
    message(paste0("found matching transcriptome:\n[ ",
                   txomeInfo$source," - ",txomeInfo$organism," - release ",txomeInfo$release," ]"))
    return(txomeInfo)
  }
  
  return(NULL)
}

# build or load a TxDb for the dataset
getTxDb <- function(txomeInfo) {
  # TODO what if there are multiple GTF files?
  stopifnot(length(txomeInfo$gtf) == 1)
  stopifnot(txomeInfo$gtf != "")
  txdbName <- basename(txomeInfo$gtf)
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, txdbName)
  if (bfccount(q) == 0) {
    # TODO: this next line already creates an entry,
    # but will need to clean up if the TxDb construction below fails
    savepath <- bfcnew(bfc, txdbName, ext=".sqlite")
    if (txomeInfo$source == "Ensembl") {
      message("building EnsDb with 'ensembldb' package")
      # TODO suppress warnings here from GTF construction?
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

# check to see if there are any missing transcripts not available
# for the rows of the tximport assay matrices. if so, give warning and subset
# (or error if all are missing)
checkAssays2Txps <- function(assays, txps) {
  assay.nms <- rownames(assays[["counts"]])
  txps.missing <- !assay.nms %in% names(txps)
  if (!all(assay.nms %in% names(txps))) {
    
    if (all(!assay.nms %in% names(txps))) {
      stop("none of the transcripts in the quantification files are in the GTF")
    } else {
      # TODO what to do here, GTF is missing some txps in FASTA for Ensembl
      warning(paste("missing some transcripts!
",sum(txps.missing), "out of", nrow(assays[["counts"]]),
"are missing from the GTF and dropped from SummarizedExperiment output"))

      # after warning, then subset
      for (nm in names(assays)) {
        assays[[nm]] <- assays[[nm]][!txps.missing,,drop=FALSE]
      }
      
    }
  }
  assays
}

makeUnrangedSE <- function(txi, coldata, metadata) {
  assays <- txi[c("counts","abundance","length")]
  # if there are inferential replicates
  if ("infReps" %in% names(txi)) {
    infReps <- rearrangeInfReps(txi)
    assays <- c(assays, infReps)
  } else if ("variance" %in% names(txi)) {
    assays <- c(assays, txi["variance"])
  }
  SummarizedExperiment(assays=assays,
                       colData=coldata,
                       metadata=metadata)
}

rearrangeInfReps <- function(txi) {
  nreps <- ncol(txi$infReps[[1]])
  stopifnot(all(sapply(txi$infReps, ncol) == nreps))
  getCols <- function(j,l) do.call(cbind, lapply(seq_along(l), function(k)  l[[k]][,j]))
  infReps <- lapply(seq_len(nreps), getCols, txi$infReps)
  names(infReps) <- paste0("infRep",seq_len(nreps))
  infReps
}
