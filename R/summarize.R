#' Summarize transcript-level quantification to gene-level
#'
#' @param se a SummarizedExperiment
#'
#' @return a SummarizedExperiment
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
  gse <- SummarizedExperiment(assays=txi.gene[c("abundance","counts","length")],
                              rowRanges=g,
                              colData=colData(se),
                              metadata=metadata(se))
  gse
}
