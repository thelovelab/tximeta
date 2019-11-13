#' Make and load linked transcriptomes ("linkedTxome")
#'
#' For now, for details please see the vignette \code{inst/script/linked.Rmd}
#' 
#' @param indexDir the local path to the Salmon index
#' @param source the source of transcriptome (e.g. "GENCODE", "Ensembl", "de-novo")
#' @param organism organism (e.g. "Homo sapiens")
#' @param release release number (e.g. "27")
#' @param genome genome (e.g. "GRCh38", or "none")
#' @param fasta location(s) for the FASTA transcript sequences
#' (of which the transcripts used to build the index is equal or a subset).
#' This can be a local path, or an HTTP or FTP URL
#' @param gtf location for the GTF/GFF file
#' (of which the transcripts used to build the index is equal or a subset).
#' This can be a local path, or an HTTP or FTP URL
#' While the \code{fasta} argument can take a vector of length greater than one
#' (more than one FASTA file containing transcripts used in indexing),
#' the \code{gtf} argument has to be a single GTF/GFF file. If transcripts
#' were added to a standard set of reference transcripts (e.g. fusion genes,
#' or pathogen transcripts), it is recommended that the tximeta user would manually
#' add these to the GTF/GFF file, and post the modified GTF/GFF publicly, such as
#' on Zenodo. This enables consistent annotation and downstream annotation
#' tasks, such as by \code{summarizeToGene}.
#' @param write logical, should a JSON file be written out
#' which documents the transcriptome checksum and metadata? (default is TRUE)
#' @param jsonFile the path to the json file for the linkedTxome
#'
#' @return nothing, the function is run for its side effects
#' 
#' @name linkedTxome
#' @rdname linkedTxome
#'
#' @examples
#'
#' # point to a Salmon quantification file with an additional artificial transcript
#' dir <- system.file("extdata/salmon_dm", package="tximportData")
#' file <- file.path(dir, "SRR1197474.plus", "quant.sf")
#' coldata <- data.frame(files=file, names="SRR1197474", sample="1",
#'                       stringsAsFactors=FALSE)
#'
#' # now point to the Salmon index itself to create a linkedTxome
#' # as the index will not match a known txome
#' indexDir <- file.path(dir, "Dm.BDGP6.22.98.plus_salmon-0.14.1")
#'
#' # point to the source FASTA and GTF:
#' fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
#'               "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz",
#'               "extra_transcript.fa.gz")
#' gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.plus.gtf.gz")
#'
#' # now create a linkedTxome, linking the Salmon index to its FASTA and GTF sources
#' makeLinkedTxome(indexDir=indexDir, source="Ensembl", organism="Drosophila melanogaster",
#'                 release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#'
#' # to clear the entire linkedTxome table
#' # (don't run unless you want to clear this table!)
#' # bfcloc <- getTximetaBFC()
#' # bfc <- BiocFileCache(bfcloc)
#' # bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
#' 
#' @export
makeLinkedTxome <- function(indexDir, source, organism, release,
                            genome, fasta, gtf, write=TRUE, jsonFile) {
  indexJson <- file.path(indexDir, "info.json")
  if (!file.exists(indexJson)) {
    indexJson <- file.path(indexDir, "header.json")
  }
  indexList <- fromJSON(indexJson)
  # Salmon's SHA-256 hash of the index is called "SeqHash" in the index JSON
  indexSeqHash <- indexList$value0$SeqHash
  # here and in the data frame where we record linkedTxome's,
  # 'index' is just the basename of the Salmon index
  index <- basename(indexDir)
  # standardize capitalization
  std.sources <- c("GENCODE","Ensembl")
  for (src in std.sources) {
    if (tolower(source) == tolower(src)) {
      source <- src
    }
  }
  # a single-row tibble for the linkedTxomeTbl
  lt <- tibble(index=index,
               source=source,
               organism=organism,
               release=release,
               genome=genome,
               fasta=list(fasta),
               gtf=gtf,
               sha256=indexSeqHash)
  stopifnot(nrow(lt) == 1)
  if (write) {
    if (missing(jsonFile)) {
      jsonFile <- paste0(indexDir,".json")
    }
    message(paste("writing linkedTxome to", jsonFile))
    # TODO be more careful about writing to a file (ask)
    write(toJSON(lt, pretty=TRUE), file=jsonFile)
  }
  stashLinkedTxome(lt)
}

#' @name linkedTxome
#' @rdname linkedTxome
#' 
#' @export
loadLinkedTxome <- function(jsonFile) {
  stashLinkedTxome(do.call(tibble, fromJSON(jsonFile)))
}

# given a single-row tibble 'lt', save this into the linkedTxomeTbl
# (the linkedTxome tibble lives in the tximeta BiocFileCache)
stashLinkedTxome <- function(lt) {
  stopifnot(is(lt, "tbl"))
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, "linkedTxomeTbl")
  if (bfccount(q) == 0) {
    message("saving linkedTxome in bfc (first time)")
    savepath <- bfcnew(bfc, "linkedTxomeTbl", ext=".rds")
    linkedTxomeTbl <- lt
    saveRDS(linkedTxomeTbl, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, "linkedTxomeTbl")
    linkedTxomeTbl <- readRDS(loadpath)
    if (lt$index %in% linkedTxomeTbl$index) {
      m <- match(lt$index, linkedTxomeTbl$index)
      stopifnot(length(m) == 1)
      if (all(mapply(identical, lt, linkedTxomeTbl[m,]))) {
        message("linkedTxome is same as already in bfc")
      } else {
        message("linkedTxome was different than one in bfc, replacing")
        linkedTxomeTbl[m,] <- lt
      }
    } else {
      message("saving linkedTxome in bfc")
      linkedTxomeTbl <- rbind(linkedTxomeTbl, lt)
    }
    saveRDS(linkedTxomeTbl, file=loadpath)
  }
  invisible()
}
