#' Add IDs to rowRanges of a SummarizedExperiment
#'
#' For now this just works with SummarizedExperiments
#' with Ensembl gene or transcript IDs. See example
#' of usage in tximeta vignette. For obtaining
#' multiple matching IDs for each row of the SummarizedExperiment
#' set \code{multiVals="list"}. See \code{select} for documentation
#' on use of \code{multiVals}.
#' 
#' @param se the SummarizedExperiment
#' @param column the name of the new ID to add (a \code{column} of the org database)
#' @param gene logical, whether the rows are genes or transcripts (default is FALSE)
#' @param ... arguments passed to \code{mapIds}
#'
#' @return a SummarizedExperiment
#'
#' @examples
#'
#' example(tximeta)
#' library(org.Dm.eg.db)	
#' se <- addIds(se, "REFSEQ", gene=FALSE)
#' 
#' @export
addIds <- function(se, column, gene=FALSE, ...) {
  missingMetadata(se, summarize=FALSE)
  # here a hack, for now this is just a prototype for ENSEMBL genes or txps
  stopifnot(metadata(se)$txomeInfo$source %in% c("Gencode","Ensembl"))
  # TODO probably should switch from package-based 'org.Hs.eg.db' to the
  # instead obtain and load OrgDb through AnnotationHub?
  orgpkg.name <- paste0("org.",
                        sub("(.).* (.).*","\\1\\2",metadata(se)$txomeInfo$organism),
                        ".eg.db")
  stopifnot(requireNamespace(orgpkg.name, quietly=TRUE))
  if (!isNamespaceLoaded(orgpkg.name)) attachNamespace(orgpkg.name)
  orgpkg <- get(orgpkg.name)
  message("mapping to new IDs using '", orgpkg.name, "' data package
if all matching IDs are desired, and '1:many mappings' are reported,
set multiVals='list' to obtain all the matching IDs")
  if (!gene & grepl("ENSG",rownames(se)[1])) {
    message("it appears the rows are gene IDs, setting 'gene' to TRUE")
    gene <- TRUE
  }
  keytype <- if (gene) "ENSEMBL" else "ENSEMBLTRANS"
  newIds <- mapIds(orgpkg, sub("\\..*","", rownames(se)),
                                     column, keytype, ...)
  mcols(se)[[column]] <- newIds
  se
}
