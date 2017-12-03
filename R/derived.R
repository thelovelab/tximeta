#' Make and load derived transcriptomes
#' 
#' @param indexDir indexDir
#' @param source source
#' @param organism organism
#' @param version version
#' @param genome genome
#' @param fasta fasta 
#' @param gtf gtf
#' @param write should a JSON file be written out which documents the transcriptome signature and metadata
#'
#' @name derivedTxome
#' @rdname derivedTxome
#' 
#' @export
makeDerivedTxome <- function(indexDir, source, organism, version, genome, fasta, gtf, write=TRUE) {
  indexList <- fromJSON(file.path(indexDir,"header.json"))
  indexSeqHash <- indexList$value0$SeqHash
  index <- basename(indexDir)
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
  
  bfc <- BiocFileCache(".")
  q <- bfcquery(bfc, "derivedTxomeDF")

  if (bfccount(q) == 0) {
    message("saving derivedTxome in bfc (first time)")
    savepath <- bfcnew(bfc, "derivedTxomeDF", ext="rds")
    derivedTxomeDF <- data.frame(dt, stringsAsFactors=FALSE)
    saveRDS(derivedTxomeDF, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, "derivedTxomeDF")
    derivedTxomeDF <- readRDS(loadpath)
    if (index %in% derivedTxomeDF$index) {
      m <- match(index, derivedTxomeDF$index)
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

# need a loadDerivedTxome() which will share a lot of code with the above function
