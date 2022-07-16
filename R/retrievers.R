#' Retrieve the TxDb or EnsDb associated with a SummarizedExperiment
#'
#' SummarizedExperiment objects returned by \code{\link{tximeta}} have
#' associated TxDb or EnsDb databases which are cached locally and used
#' to perform various metadata related tasks. This helper function
#' retrieves the database itself for the user to perform any additional
#' operations.
#'
#' @param se the SummarizedExperiment
#'
#' @return a database object
#'
#' @examples
#'
#' example(tximeta)
#' edb <- retrieveDb(se)
#' 
#' @export
retrieveDb <- function(se) {
  missingMetadata(se, summarize=FALSE)
  txomeInfo <- metadata(se)$txomeInfo
  getTxDb(txomeInfo)
}

#' Retrieve the cDNA transcript sequence for a SummarizedExperiment
#'
#' This helper function retrieves the cDNA sequence of
#' the transcripts used for expression quantification.
#' This function either downloads or loads the transcript
#' sequence from cache, it does not re-order or check against
#' the rows of the SummarizedExperiment (which could be
#' already summarized to genes for example).
#'
#' @param se the SummarizedExperiment
#' @param quiet logical, suppress messages
#'
#' @return a DNAStringSet object
#'
#' @examples
#'
#' \dontrun{
#' # this example is not run because it requires access to Ensembl ftp
#' example(tximeta)
#' cdna <- retrieveCDNA(se)
#' }
#' 
#' @export
retrieveCDNA <- function(se, quiet=FALSE) {
  if (!is(se, "SummarizedExperiment")) {
    txomeInfo <- se
  } else {
    missingMetadata(se, summarize=FALSE)
    txomeInfo <- metadata(se)$txomeInfo
  }
  cdna.id <- paste0("cdna-",substr(txomeInfo$sha256,1,32))
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, cdna.id)
  if (bfccount(q) == 1) {
    loadpath <- bfcrpath(bfc, cdna.id)
    cdna <- readRDS(loadpath)
  } else { 
    if (is.list(txomeInfo$fasta)) {
      cdna <- list()
      for (i in seq_along(txomeInfo$fasta[[1]])) {
        cdna[[i]] <- readDNAStringSet(txomeInfo$fasta[[1]][i])
      }
      cdna <- do.call(c, cdna)
    } else {
      cdna <- readDNAStringSet(txomeInfo$fasta)
    }
    # save it to BFC so as not to repeat the above operation
    savepath <- bfcnew(bfc, cdna.id, ext=".rds")
    saveRDS(cdna, file=savepath)
  }
  if (!quiet) message("retrieving cDNA, note that sequence is not ordered or matched to the rows of 'se'")
  cdna
}
