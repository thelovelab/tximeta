#' Add IDs to rowRanges of a SummarizedExperiment
#'
#' For now this just works with SummarizedExperiments with Ensembl gene or transcript IDs
#' 
#' @param se the SummarizedExperiment
#' @param column the name of the new ID to add
#'
#' @return a SummarizedExperiment
#'
#' @export
addIds <- function(se, column) {
  # here a hack, for now this is just a prototype for ENSEMBL genes or txps
  stopifnot(metadata(se)$txomeInfo$source %in% c("Gencode","Ensembl"))
  isGene <- substr(rownames(se)[1],1,4) == "ENSG"
  orgpkg.name <- sub(" ",".",metadata(se)$txomeInfo$organism)
  stopifnot(requireNamespace(orgpkg.name, quietly=TRUE))
  orgpkg <- get(orgpkg.name)
  message(paste0("mapping to new IDs using '", orgpkg.name, "' data package"))
  keytype <- if (isGene) "ENSEMBL" else "ENSEMBLTRANS"
  suppressMessages({newIds <- mapIds(orgpkg, sub("\\..*","", rownames(se)),
                                     column, keytype)})
  mcols(se)[[column]] <- newIds
  se
}
