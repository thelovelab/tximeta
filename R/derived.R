#' Make and load derived transcriptomes ("derivedTxome")
#'
#' For now, for details please see the vignette \code{inst/script/derived.Rmd}
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
#' @param jsonFile the path to the json file for the derivedTxome
#'
#' @name derivedTxome
#' @rdname derivedTxome
#' 
#' @export
makeDerivedTxome <- function(indexDir, source, organism, version,
                             genome, fasta, gtf, write=TRUE, jsonFile) {
  indexList <- fromJSON(file.path(indexDir,"header.json"))
  indexSeqHash <- indexList$value0$SeqHash
  # here and in the data frame where we record derivedTxome's,
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
    message(paste("writing derivedTxome to", jsonFile))
    # TODO be more careful about writing to a file (ask)
    write(toJSON(dt, pretty=TRUE), file=jsonFile)
  }
  stashDerivedTxome(dt)
}

#' @name derivedTxome
#' @rdname derivedTxome
#' 
#' @export
loadDerivedTxome <- function(jsonFile) {
  stashDerivedTxome(do.call(tibble, fromJSON(jsonFile)))
}

stashDerivedTxome <- function(dt) {
  stopifnot(is(dt, "tbl"))
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, "derivedTxomeDF")
  if (bfccount(q) == 0) {
    message("saving derivedTxome in bfc (first time)")
    savepath <- bfcnew(bfc, "derivedTxomeDF", ext="rds")
    derivedTxomeDF <- dt
    saveRDS(derivedTxomeDF, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, "derivedTxomeDF")
    derivedTxomeDF <- readRDS(loadpath)
    if (dt$index %in% derivedTxomeDF$index) {
      m <- match(dt$index, derivedTxomeDF$index)
      stopifnot(length(m) == 1)
      if (all(mapply(all.equal, dt, derivedTxomeDF[m,]))) {
        message("derivedTxome is same as already in bfc")
      } else {
        message("derivedTxome was different than one in bfc, replacing")
        derivedTxomeDF[m,] <- dt
      }
    } else {
      message("saving derivedTxome in bfc")
      derivedTxomeDF <- rbind(derivedTxomeDF, dt)
    }
  }
  invisible()
}
