#' Add IDs to rowRanges of a SummarizedExperiment
#'
#' For now this function just works with SummarizedExperiments
#' with Ensembl gene or transcript IDs. See example
#' of usage in tximeta vignette. For obtaining
#' multiple matching IDs for each row of the SummarizedExperiment
#' set \code{multiVals="list"}. See \code{select} for documentation
#' on use of \code{multiVals}.
#' 
#' @param se the SummarizedExperiment
#' @param column the name of the new ID to add (a \code{column} of the org package
#' database or of the TxDb/EnsDb is \code{fromDb=TRUE})
#' @param fromDb logical, whether to use the TxDb/EnsDb that is associated
#' with \code{se}. Default is FALSE, and an org package is used.
#' Currently only implemented for transcript level (gene=FALSE).
#' Column names can be viewed with \code{columns(retrieveDb(se))}
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
addIds <- function(se, column, fromDb=FALSE, gene=FALSE, ...) {
  missingMetadata(se, summarize=FALSE)
  keys <- rownames(se)
  if (!fromDb) {
    # works for GENCODE or Ensembl, or linkedTxomes with local GTF supplying the DB
    allowedSourceNames <- c("GENCODE","Ensembl","LocalGENCODE","LocalEnsembl")
    stopifnot(metadata(se)$txomeInfo$source %in% allowedSourceNames)
    # TODO probably should switch from package-based 'org.Hs.eg.db'
    # and instead obtain and load OrgDb through AnnotationHub
    dbname <- paste0("org.",
                     sub("(.).* (.).*","\\1\\2",metadata(se)$txomeInfo$organism),
                     ".eg.db")
    stopifnot(requireNamespace(dbname, quietly=TRUE))
    if (!isNamespaceLoaded(dbname)) attachNamespace(dbname)
    db <- get(dbname)
    if (!gene & grepl("ENSG|ENSMUSG",rownames(se)[1])) {
      message("it appears the rows are gene IDs, setting 'gene' to TRUE")
      gene <- TRUE
    }
    if (gene & grepl("ENST|ENSMUST",rownames(se)[1])) {
      message("gene=TRUE and rows are transcripts: using 'gene_id' column to map IDs")
      stopifnot("gene_id" %in% names(mcols(se)))
      keys <- mcols(se)$gene_id
    }
    keytype <- if (gene) "ENSEMBL" else "ENSEMBLTRANS"
    keys <- sub("\\..*","", keys)
  } else {
    if (gene) stop("fromDb=TRUE implemented only for gene=FALSE")
    db <- getTxDb(metadata(se)$txomeInfo)
    dbname <- class(db)
    keytype <- "TXNAME"
  }
  message(paste0("mapping to new IDs using ", dbname, "
if all matching IDs are desired, and '1:many mappings' are reported,
set multiVals='list' to obtain all the matching IDs"))
  newIds <- mapIds(db, keys, column, keytype, ...)
  mcols(se)[[column]] <- newIds
  se
}
