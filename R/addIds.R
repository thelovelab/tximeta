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
#' @param gene logical, whether to map by genes or transcripts (default is FALSE).
#' if rows are genes, and easily detected as such (ENSG or ENSMUSG), it will
#' automatically switch to TRUE. if rows are transcripts and \code{gene=TRUE},
#' then it will try to use a \code{gene_id} column to map IDs to \code{column}
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
  if (!gene & grepl("ENSG|ENSMUSG",rownames(se)[1])) {
    message("it appears the rows are gene IDs, setting 'gene' to TRUE")
    gene <- TRUE
  }
  keys <- rownames(se)
  if (gene & grepl("ENST|ENSMUST",rownames(se)[1])) {
    message("gene=TRUE and rows are transcripts: using 'gene_id' column to map IDs")
    stopifnot("gene_id" %in% names(mcols(se)))
    keys <- mcols(se)$gene_id
  }
  keytype <- if (gene) "ENSEMBL" else "ENSEMBLTRANS"
  newIds <- mapIds(orgpkg, sub("\\..*","", keys),
                   column, keytype, ...)
  mcols(se)[[column]] <- newIds
  se
}
