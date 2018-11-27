summarizeToGene.SummarizedExperiment <- function(object, varReduce=FALSE, ...) {

  missingMetadata(object, summarize=TRUE)
    
  txdb <- getTxDb(metadata(object)$txomeInfo)
  message("obtaining transcript-to-gene mapping from TxDb")

  # TODO fix this next line of code for EnsDb
  suppressMessages({tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")})

  # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?
  suppressWarnings({ g <- genes(txdb) })

  txi <- list(
    abundance=assays(object)[["abundance"]],
    counts=assays(object)[["counts"]],
    length=assays(object)[["length"]],
    countsFromAbundance="no"
  )
  if ("infRep1" %in% assayNames(object)) {
    infReps <- list(assays(object)[grep("infRep", assayNames(object))])
    if (varReduce) {
      # split from per replicate list into per sample list
      # (this is what tximport expects for varReduce=TRUE)
      infReps <- list(splitInfReps(infReps[[1]]))
    }
    txi <- c(txi, infReps=infReps)
  } else {
    if (varReduce) stop("cannot calculate inferential variance without inferential replicates")
  }

  txi.gene <- summarizeToGene(object=txi, tx2gene=tx2gene, varReduce=varReduce, ...)

  # put 'counts' in front to facilitate DESeqDataSet construction
  assays <- txi.gene[c("counts","abundance","length")]
  if (varReduce) {
    assays <- c(assays, txi.gene["variance"])
  } else if ("infRep1" %in% assayNames(object)) {
    infReps <- txi.gene$infReps
    assays <- c(assays, infReps)
  }

  # here do the same check/subset but with gene-level tximport assay matrices and gene ranges
  assays <- checkAssays2Txps(assays, g)

  # TODO give a warning here if there are genes in TxDb not in Salmon index?
  g <- g[rownames(assays[["counts"]])]

  metadata <- metadata(object)
  # stash countsFromAbundance value
  metadata$countsFromAbundance <- txi.gene$countsFromAbundance
  
  gse <- SummarizedExperiment(assays=assays,
                              rowRanges=g,
                              colData=colData(object),
                              metadata=metadata)  
  gse
}

#' Summarize estimated quantitites to gene-level
#'
#' Summarizes abundances, counts, lengths, (and inferential
#' replicates or variance) from transcript- to gene-level.
#' This function operates on SummarizedExperiment objects, and
#' will automatically access the relevant TxDb (by either finding it
#' in the BiocFileCache or by building it from an ftp location).
#' #' This function uses the tximport package to perform summarization,
#' where a method is defined that works on simple lists.
#'
#' @param object a SummarizedExperiment produced by \code{tximeta}
#' @param varReduce whether to reduce per-sample inferential replicates
#' information into a matrix of sample variances \code{variance} (default FALSE)
#' @param ... arguments passed to \code{tximport}
#'
#' @return a SummarizedExperiment with summarized quantifications
#'
#' @rdname summarizeToGene
#' @docType methods
#' @aliases summarizeToGene,SummarizedExperiment-method
#' 
#' @examples
#'
#' example(tximeta)
#' gse <- summarizeToGene(se)
#'
#' @importFrom tximport summarizeToGene
#' 
#' @export
setMethod("summarizeToGene", signature(object="SummarizedExperiment"),
          summarizeToGene.SummarizedExperiment)

missingMetadata <- function(se, summarize=TRUE) {
  msg <- "use of this function requires transcriptome metadata which is missing.
  use a linkedTxome to provide the missing metadata and rerun tximeta()"
  if (summarize) {
    msg <- paste0(msg, "
  or use tx2gene, txOut=FALSE (and skipMeta=TRUE if Salmon/Sailfish)")
  }
  if (is.null(metadata(se)$txomeInfo)) stop(msg)
}
