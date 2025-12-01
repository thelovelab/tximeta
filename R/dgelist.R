#' Make a DGEList from tximeta output
#'
#' A simple wrapper function for constructing a DGEList for use with edgeR.
#' See vignette for an example. Requires installation of the edgeR
#' package from Bioconductor.
#'
#' @param se a SummarizedExperiment produced by tximeta
#'
#' @return a DGEList
#'
#' @export
makeDGEList <- function(se) {
  if (!requireNamespace("edgeR", quietly=TRUE)) {
    stop("this function requires the edgeR package is installed")
  }
  cts <- assays(se)[["counts"]]
  normMat <- assays(se)[["length"]]
  normMat <- normMat / exp(rowMeans(log(normMat)))
  o <- log(edgeR::calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  y <- edgeR::DGEList(cts, samples=as.data.frame(colData(se)),
                      genes=as.data.frame(rowData(se)))
  y <- edgeR::scaleOffset(y, t(t(log(normMat)) + o))
  y
}
