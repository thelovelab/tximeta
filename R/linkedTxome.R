#' Make and load linked transcriptomes ("linkedTxome")
#'
#' For now, for details please see the vignette \code{inst/script/linked.Rmd}
#' 
#' @param indexDir the path to the Salmon or Sailfish index
#' @param source the source of transcriptome (e.g. "Gencode" or "Ensembl")
#' @param organism organism (e.g. "Homo sapiens")
#' @param version version (e.g. "27")
#' @param genome genome (e.g. "GRCh38")
#' @param fasta FTP location for the FASTA sequence (of which the index is a subset)
#' @param gtf FTP location for the GTF file (of which the index is a subset)
#' @param write should a JSON file be written out
#' which documents the transcriptome signature and metadata? (default is TRUE)
#' @param jsonFile the path to the json file for the linkedTxome
#'
#' @name linkedTxome
#' @rdname linkedTxome
#' 
#' @export
makeLinkedTxome <- function(indexDir, source, organism, version,
                             genome, fasta, gtf, write=TRUE, jsonFile) {
  indexList <- fromJSON(file.path(indexDir,"header.json"))
  indexSeqHash <- indexList$value0$SeqHash
  # here and in the data frame where we record linkedTxome's,
  # 'index' is just the basename of the Salmon index
  index <- basename(indexDir)
  # standardize capitalization
  std.sources <- c("Gencode","Ensembl")
  for (src in std.sources) {
    if (tolower(source) == tolower(src)) {
      source <- src
    }
  }
  dt <- tibble(index=index,
               index_seq_hash=indexSeqHash,
               source=source,
               organism=organism,
               version=version,
               genome=genome,
               fasta=list(fasta),
               gtf=gtf)
  stopifnot(nrow(dt) == 1)
  if (write) {
    if (missing(jsonFile)) {
      jsonFile <- paste0(indexDir,".json")
    }
    message(paste("writing linkedTxome to", jsonFile))
    # TODO be more careful about writing to a file (ask)
    write(toJSON(dt, pretty=TRUE), file=jsonFile)
  }
  stashLinkedTxome(dt)
}

#' @name linkedTxome
#' @rdname linkedTxome
#' 
#' @export
loadLinkedTxome <- function(jsonFile) {
  stashLinkedTxome(do.call(tibble, fromJSON(jsonFile)))
}

stashLinkedTxome <- function(dt) {
  stopifnot(is(dt, "tbl"))
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, "linkedTxomeTbl")
  if (bfccount(q) == 0) {
    message("saving linkedTxome in bfc (first time)")
    savepath <- bfcnew(bfc, "linkedTxomeTbl", ext="rds")
    linkedTxomeTbl <- dt
    saveRDS(linkedTxomeTbl, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, "linkedTxomeTbl")
    linkedTxomeTbl <- readRDS(loadpath)
    if (dt$index %in% linkedTxomeTbl$index) {
      m <- match(dt$index, linkedTxomeTbl$index)
      stopifnot(length(m) == 1)
      if (all(mapply(all.equal, dt, linkedTxomeTbl[m,]))) {
        message("linkedTxome is same as already in bfc")
      } else {
        message("linkedTxome was different than one in bfc, replacing")
        linkedTxomeTbl[m,] <- dt
      }
    } else {
      message("saving linkedTxome in bfc")
      linkedTxomeTbl <- rbind(linkedTxomeTbl, dt)
    }
  }
  invisible()
}
