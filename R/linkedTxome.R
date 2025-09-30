#' Make and load linked transcriptomes (linked GTF and FASTA)
#'
#' `makeLinkedTxome()` reads the digest associated with a salmon
#' index at `indexDir`, and persistently links it to metadata
#' (alternatively the `digest` string itself and an 
#' `indexName` can be provided).
#' Linked metadata includes key information
#' about the transcriptome, including the `source`, `organism`,
#' `release`, and `genome` (these are custom character strings),
#' as well as the locations (e.g. local, HTTP, or FTP) for one or more `fasta`
#' files and one `gtf` file. `loadLinkedTxome()` loads this
#' information from a JSON file. See _Details_.
#'
#' `makeLinkedTxome()` links the information about the transcriptome
#' used for quantification in two ways:
#' 1) the function will store a record in tximeta's cache such that
#' future import of quantification data will automatically access and
#' parse the GTF as if the transcriptome were one of those automatically
#' detected by tximeta. Then all features of tximeta (e.g. summarization
#' to gene, programmatic adding of IDs or metadata) will be available;
#' 2) it will by default write out a JSON file
#' that can be shared, or posted online, and which can be read by
#' `loadLinkedTxome()` which will store the information in tximeta's
#' cache. This should make the full quantification-import pipeline
#' computationally reproducible / auditable even for transcriptomes
#' which differ from those provided by references (GENCODE, Ensembl,
#' RefSeq).
#' 
#' For further details please see the "Linked transcriptomes"
#' section of the tximeta vignette.
#' 
#' This function can be used in combination with `tximixInspectDigests()`
#' for use with `tximix()`-imported oarfish data, when multiple
#' reference transcript sets have been indexed.
#' 
#' @param digest the full digest as character string,
#' (this or `indexDir` is required, only one should be specified)
#' @param indexName a name for the `index` when storing the linkedTxome,
#' required if providing the `digest` string, suggest using the
#' basename of the FASTA file and the software used, 
#' e.g. "gencode.vXX_salmon-0.XX.Y"
#' @param indexDir the local path to the salmon index 
#' (this or `digest` is required, only one should be specified)
#' @param source the source of transcriptome (e.g. "de-novo").
#' Note: if you specify "GENCODE" or "Ensembl", this will trigger
#' behavior by tximeta that may not be desired: e.g. attempts to
#' download canonical transcriptome data from AnnotationHub
#' (unless useHub=FALSE when running tximeta) and parsing of
#' Ensembl GTF using ensembldb (which may fail if the GTF file
#' has been modified). For transcriptomes that are defined by
#' local GTF files, it is recommended to use the terms "LocalGENCODE"
#' or "LocalEnsembl". Setting "LocalEnsembl" will also strip
#' version numbers from the FASTA transcript IDs to enable matching
#' with the Ensembl GTF.
#' @param organism organism (e.g. "Homo sapiens")
#' @param release release number (e.g. "27")
#' @param genome genome (e.g. "GRCh38", or "none")
#' @param fasta location(s) for the FASTA transcript sequences
#' (of which the transcripts used to build the index is equal or a subset).
#' This can be a local path, or an HTTP or FTP URL
#' @param gtf location for the GTF/GFF file
#' (of which the transcripts used to build the index is equal or a subset).
#' This can be a local path, or an HTTP or FTP URL
#' While the `fasta` argument can take a vector of length greater than one
#' (more than one FASTA file containing transcripts used in indexing),
#' the `gtf` argument has to be a single GTF/GFF file.
#' This can also be a serialized GRanges object (location of a .rds file)
#' imported with rtracklayer.
#' If transcripts were added to a standard set of reference transcripts (e.g. fusion genes,
#' or pathogen transcripts), it is recommended that the tximeta user would manually
#' add these to the GTF/GFF file, and post the modified GTF/GFF publicly, such as
#' on Zenodo. This enables consistent annotation and downstream annotation
#' tasks, such as by
#' [`summarizeToGene()`][summarizeToGene,SummarizedExperiment-method].
#' @param write logical, should a JSON file be written out
#' which documents the transcriptome digest and metadata? (default is TRUE)
#' @param jsonFile the path to the json file for the linkedTxome
#'
#' @return nothing, the function is run for its side effects
#' 
#' @name linkedTxome
#' @rdname linkedTxome
#'
#' @examples
#'
#' # point to a salmon quantification file with an additional artificial transcript
#' dir <- system.file("extdata/salmon_dm", package="tximportData")
#' file <- file.path(dir, "SRR1197474.plus", "quant.sf")
#' coldata <- data.frame(files=file, names="SRR1197474", sample="1",
#'                       stringsAsFactors=FALSE)
#'
#' # now point to the salmon index itself to create a linkedTxome
#' # as the index will not match a known txome
#' indexDir <- file.path(dir, "Dm.BDGP6.22.98.plus_salmon-0.14.1")
#'
#' # point to the source FASTA and GTF:
#' fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
#'               "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz",
#'               "extra_transcript.fa.gz")
#' gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.plus.gtf.gz")
#'
#' # now create a linkedTxome, linking the salmon index to its FASTA and GTF sources
#' makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster",
#'                 release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#'
#' # to clear the entire linkedTxome table
#' # (don't run unless you want to clear this table!)
#' # bfcloc <- getTximetaBFC()
#' # bfc <- BiocFileCache(bfcloc)
#' # bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
#' 
#' @export
makeLinkedTxome <- function(
  digest=NULL,
  indexName,
  indexDir=NULL,
  source, organism, release,
  genome, fasta, gtf, write=TRUE, jsonFile
) {
  
  # only one or the other is specified
  stopifnot(xor(is.null(digest), is.null(indexDir)))

  if (!is.null(indexDir)) {
    # `indexDir` was specified
    message(paste0("reading digest from indexDir: ",indexDir))
    indexJson <- file.path(indexDir, "info.json")
    # backup spot for information...
    if (!file.exists(indexJson)) {
      indexJson <- file.path(indexDir, "header.json")
    }
    indexList <- fromJSON(indexJson)
    # salmon's SHA-256 hash of the index is called "SeqHash" in the index JSON
    # Pre-salmon 1.0.0 the header.json file has a "value0" sublist, 
    # from salmon 1.0.0 the info.json file doesn't
    if ("value0" %in% names(indexList)) {
      digest <- indexList$value0$SeqHash
    } else {
      digest <- indexList$SeqHash
    }
    # here and in the data frame where we record linkedTxome's,
    # 'index' is just the basename of the salmon index
    index <- basename(indexDir)
  } else {
    # `digest` was specified, so use the indexName provided
    stopifnot(!missing(indexName))
    message(paste0("linking file-based metadata to digest: ",substr(1,6,digest),"..."))
    index <- indexName
  }

  # standardize capitalization
  std.sources <- c("GENCODE","Ensembl")
  for (src in std.sources) {
    if (tolower(source) == tolower(src)) {
      source <- src
    }
  }
  if (source %in% std.sources) {
    if (source == "Ensembl") {
      message("NOTE: linkedTxome with source='Ensembl', ensembldb will be used to parse GTF.
this may produce errors if the GTF is not from Ensembl, or has been modified.
set useHub=FALSE in tximeta to avoid download of reference txome from AnnotationHub.
alternatively use a different string for source argument")
    } else {
      message("NOTE: linkedTxome with source='GENCODE', set useHub=FALSE in tximeta
to avoid download of reference txome from AnnotationHub.
alternatively use a different string for source argument")
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
               sha256=digest)
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

### linkedTxpData -- a lightweight alternative ###

#' Make linked transcript data (linked GRanges)
#' 
#' @param digest character string of the full digest of the 
#' reference transcripts, see `tximixInspectDigests()` with `fullDigest=TRUE`
#' @param txpData _GRanges_ providing information about ranges 
#' representing the transcript sequences linked to `digest`
#' 
#' @return nothing, the function is run for its side effects
#' 
#' @name linkedTxpData
#' @rdname linkedTxpData
#' 
#' @export
makeLinkedTxpData <- function(digest, txpData) {
    message(paste0("linking user-provided metadata to digest: ",substr(digest,1,6),"..."))
}
