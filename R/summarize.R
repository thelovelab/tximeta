#' Summarize transcript-level quantifications to gene-level
#'
#' This function uses the tximport package to summarize counts,
#' abundances and effective lengths from transcript-level to gene-level.
#' It automatically will access the relevant TxDb (by either finding it
#' in the BiocFileCache or by building it from an ftp location).
#'
#' @param se a SummarizedExperiment produced by \code{tximeta}
#' @param varReduce whether to reduce per-sample inferential replicates
#' information into a matrix of sample variances \code{variance} (default FALSE)
#' @param ... arguments passed to \code{tximport}
#'
#' @return a SummarizedExperiment with summarized quantifications
#'
#' @examples
#'
#' example(tximeta)
#' gse <- summarizeToGene(se)
#'
#' @export
summarizeToGene <- function(se, varReduce=FALSE, ...) {

  # TODO make `summarizeToGene` a generic in tximport

  txdb <- getTxDb(metadata(se)$txomeInfo)
  message("obtaining transcript-to-gene mapping from TxDb")

  # TODO fix this next line of code for EnsDb
  suppressMessages({tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")})

  # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?
  suppressWarnings({ g <- genes(txdb) })

  txi <- list(
    abundance=assays(se)[["abundance"]],
    counts=assays(se)[["counts"]],
    length=assays(se)[["length"]],
    countsFromAbundance="no"
  )
  if ("infRep1" %in% assayNames(se)) {
    infReps <- list(assays(se)[grep("infRep", assayNames(se))])
    if (varReduce) {
      # split from per replicate list into per sample list
      # (this is what tximport expects for varReduce=TRUE)
      infReps <- list(splitInfReps(infReps[[1]]))
    }
    txi <- c(txi, infReps=infReps)
  } else {
    if (varReduce) stop("cannot calculate inferential variance without inferential replicates")
  }

  txi.gene <- tximport::summarizeToGene(txi, tx2gene, varReduce=varReduce, ...)

  # put 'counts' in front to facilitate DESeqDataSet construction
  assays <- txi.gene[c("counts","abundance","length")]
  if (varReduce) {
    assays <- c(assays, txi.gene["variance"])
  } else if ("infRep1" %in% assayNames(se)) {
    infReps <- txi.gene$infReps
    assays <- c(assays, infReps)
  }

  # here do the same check/subset but with gene-level tximport assay matrices and gene ranges
  assays <- checkAssays2Txps(assays, g)

  # TODO give a warning here if there are genes in TxDb not in Salmon index?
  g <- g[rownames(assays[["counts"]])]

  metadata <- metadata(se)
  # stash countsFromAbundance value
  metadata$countsFromAbundance <- txi.gene$countsFromAbundance
  
  gse <- SummarizedExperiment(assays=assays,
                              rowRanges=g,
                              colData=colData(se),
                              metadata=metadata)  
  gse
}
