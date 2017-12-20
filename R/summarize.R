#' Summarize transcript-level quantifications to gene-level
#'
#' This function uses the tximport package to summarize counts,
#' abundances and effective lengths from transcript-level to gene-level.
#' It automatically will access the relevant TxDb (by either finding it
#' in the BiocFileCache or by building it from an ftp location).
#'
#' @param se a SummarizedExperiment produced by \code{tximeta}
#'
#' @return a SummarizedExperiment with summarized quantifications
#'
#' @export
summarizeToGene <- function(se) {
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
  txi.gene <- tximport::summarizeToGene(txi, tx2gene)
  g <- genes(txdb)
  stopifnot(all(rownames(txi.gene$counts) == names(g)))
  gse <- SummarizedExperiment(assays=txi.gene[c("counts","abundance","length")],
                              rowRanges=g,
                              colData=colData(se),
                              metadata=metadata(se))
  gse
}
