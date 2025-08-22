#' Import transcript quantification across sets of transcripts
#'
#' tximix allows separation of annotated and novel transcripts
#' during transcript quantification import and addition of metadata.
#' This function supports quantification in the case of mixed
#' transcriptomes, including additions from de novo assembly,
#' as well as transgenes and spike-ins.
#'
#' @param coldata data.frame with columns \code{files} and \code{names}
#' as in \code{\link{tximeta}}
#' @param type what quantifier was used (see \code{\link{tximport}})
#' @param ... passed to tximport
#'
#' @return an unranged SummarizedExperiment
#'
#' @export
tximix <- function(coldata, type="oarfish", ...) {
  stopifnot(type == "oarfish")

  # tximeta metadata
  metadata <- makeMetadata(type)

  txi <- tximport(files, type=type, txOut=TRUE, ...)
  metadata$countsFromAbundance <- txi$countsFromAbundance

  se <- makeUnrangedSE(txi, coldata, metadata)

  return(se)
}
