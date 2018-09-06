#' Summarize transcript-level quantifications to gene-level
#'
#' This function uses the tximport package to summarize counts,
#' abundances and effective lengths from transcript-level to gene-level.
#' It automatically will access the relevant TxDb (by either finding it
#' in the BiocFileCache or by building it from an ftp location).
#'
#' @param se a SummarizedExperiment produced by \code{tximeta}
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
summarizeToGene <- function(se, ...) {

  # TODO make `summarizeToGene` a generic in tximport

  txdb <- getTxDb(metadata(se)$txomeInfo)
  message("obtaining transcript-to-gene mapping from TxDb")
  # TODO fix this next line of code for EnsDb
  suppressMessages({tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")})
  txi <- list(
    abundance=assays(se)$abundance,
    counts=assays(se)$counts,
    length=assays(se)$length,
    countsFromAbundance="no"
  )
  txi.gene <- tximport::summarizeToGene(txi, tx2gene, ...)

  # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?
  suppressWarnings({ g <- genes(txdb) })

  # here do the same check/subset but with gene-level txi matrices and gene ranges
  txi.gene <- checkTxi2Txps(txi.gene, g)

  # TODO give a warning here if there are genes in TxDb not in Salmon index?
  g <- g[rownames(txi.gene$counts)]

  metadata <- metadata(se)
  # stash countsFromAbundance value
  metadata$countsFromAbundance <- txi.gene$countsFromAbundance
  
  
  # put 'counts' in front to facilitate DESeqDataSet construction
  gse <- SummarizedExperiment(assays=txi.gene[c("counts","abundance","length")],
                              rowRanges=g,
                              colData=colData(se),
                              metadata=metadata)
  
  gse
}
