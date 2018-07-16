#' Add IDs to rowRanges of a SummarizedExperiment
#'
#' For now this just works with SummarizedExperiments
#' with Ensembl gene or transcript IDs.
#' 
#' @param se the SummarizedExperiment
#' @param column the name of the new ID to add (a \code{column} of the org database)
#' @param gene logical, whether the rows are genes or transcripts (default is FALSE)
#'
#' @return a SummarizedExperiment
#'
#' @examples
#'
#' example(tximeta)
#' library(org.Dm.eg.db)
#' # normally:
#' # se <- addIds(se, "REFSEQ", gene=FALSE)
#' 
#' # for speed we just subset to the first 10 rows:
#' se.sub <- addIds(se[1:10,], "REFSEQ", gene=FALSE)
#' 
#' @export
addIds <- function(se, column, gene=FALSE) {
  # here a hack, for now this is just a prototype for ENSEMBL genes or txps
  stopifnot(metadata(se)$txomeInfo$source %in% c("Gencode","Ensembl"))
  orgpkg.name <- paste0("org.",
                        sub("(.).* (.).*","\\1\\2",metadata(se)$txomeInfo$organism),
                        ".eg.db")
  stopifnot(requireNamespace(orgpkg.name, quietly=TRUE))
  orgpkg <- get(orgpkg.name)
  message(paste0("mapping to new IDs using '", orgpkg.name, "' data package"))
  keytype <- if (gene) "ENSEMBL" else "ENSEMBLTRANS"
  suppressMessages({newIds <- mapIds(orgpkg, sub("\\..*","", rownames(se)),
                                     column, keytype)})
  mcols(se)[[column]] <- newIds
  se
}
