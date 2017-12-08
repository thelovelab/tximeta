#' Make and load derived transcriptomes (derivedTxome)
#' 
#' @param indexDir indexDir
#' @param source source of transcriptome FASTA and GTF
#' @param organism organism
#' @param version version
#' @param genome genome
#' @param fasta fasta 
#' @param gtf gtf
#' @param write should a JSON file be written out which documents the transcriptome signature and metadata
#' @param jsonfile the path to the json file for the derivedTxome
#'
#' @name derivedTxome
#' @rdname derivedTxome
#' 
#' @export
makeDerivedTxome <- function(indexDir, source, organism, version, genome, fasta, gtf, write=TRUE) {
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
  # TODO this should be made to work for multiple fasta and GTF sources...
  dt <- list(index=index,
             index_seq_hash=indexSeqHash,
             source=source,
             organism=organism,
             version=version,
             genome=genome,
             fasta=fasta,
             gtf=gtf)
  if (write) {
    filename <- paste0(indexDir,".json")
    message(paste("writing derivedTxome to", filename))
    # TODO be more careful about writing to a file (ask)
    write(toJSON(dt, pretty=TRUE), file=filename)
  }
  stashDerivedTxome(dt)
}

#' @name derivedTxome
#' @rdname derivedTxome
#' 
#' @export
loadDerivedTxome <- function(jsonfile) {
  stashDerivedTxome(fromJSON(jsonfile))
}

stashDerivedTxome <- function(dt) {
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, "derivedTxomeDF")
  if (bfccount(q) == 0) {
    message("saving derivedTxome in bfc (first time)")
    savepath <- bfcnew(bfc, "derivedTxomeDF", ext="rds")
    derivedTxomeDF <- data.frame(dt, stringsAsFactors=FALSE)
    saveRDS(derivedTxomeDF, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, "derivedTxomeDF")
    derivedTxomeDF <- readRDS(loadpath)
    if (dt$index %in% derivedTxomeDF$index) {
      m <- match(dt$index, derivedTxomeDF$index)
      if (all(data.frame(dt) == derivedTxomeDF[m,,drop=FALSE])) {
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
