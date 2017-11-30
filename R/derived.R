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
  # TODO this should be made to work for multiple fasta and GTF sources...
  dt <- list(index_seq_hash=indexSeqHash,
             source=source,
             organism=organism,
             version=version,
             genome=genome,
             fasta=fasta,
             gtf=gtf)
  if (write) {
    filename <- paste0(indexDir,".json")
    message(paste("writing derivedTxome to", filename))
    write(toJSON(dt, pretty=TRUE), file=filename)
  }
  dt
}
