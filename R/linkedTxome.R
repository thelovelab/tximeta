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
#' This function can be used in combination with `inspectDigests()`
#' and oarfish data from `importData()`, when multiple
#' reference transcript sets have been indexed. See also
#' `makeLinkedTxpData()`.
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
#' @param prefer whether to prefer loading data from linkedTxome: `txome`,
#' linkedTxpData: `txpdata`, or the `annotated` version
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
  digest = NULL,
  indexName,
  indexDir = NULL,
  source,
  organism,
  release,
  genome,
  prefer = c("txome", "txpdata", "annotated"),
  fasta,
  gtf,
  write = TRUE,
  jsonFile
) {
  # only one or the other is specified
  stopifnot(xor(is.null(digest), is.null(indexDir)))

  if (!is.null(indexDir)) {
    # `indexDir` was specified
    message(paste0("reading digest from indexDir: ", indexDir))
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
    message(paste0(
      "linking file-based metadata to digest: ",
      substr(digest, 1, 6),
      "..."
    ))
    index <- indexName
  }

  std_sources <- c("GENCODE", "Ensembl")
  source <- standardizeCapitalization(source, std_sources)

  if (source %in% std_sources) {
    if (source == "Ensembl") {
      message(
        "NOTE: linkedTxome with source='Ensembl', ensembldb will be used to parse GTF.
this may produce errors if the GTF is not from Ensembl, or has been modified.
set useHub=FALSE in tximeta to avoid download of reference txome from AnnotationHub.
alternatively use a different string for source argument, e.g. LocalEnsembl"
      )
    } else {
      message(
        "NOTE: linkedTxome with source='GENCODE', set useHub=FALSE in tximeta
to avoid download of reference txome from AnnotationHub.
alternatively use a different string for source argument, e.g. LocalGENCODE"
      )
    }
  }
  # a single-row tibble for the linkedTxomeTbl
  # matches conent below in updateLinkedThingTbl()
  lt <- tibble(
    index = index,
    source = source,
    organism = organism,
    release = release,
    genome = genome,
    fasta = list(fasta),
    gtf = gtf,
    sha256 = digest
  )
  stopifnot(nrow(lt) == 1)
  if (write) {
    if (missing(jsonFile)) {
      jsonFile <- paste0(indexDir, ".json")
    }
    message(paste("writing linkedTxome to", jsonFile))
    # TODO be more careful about writing to a file (ask)
    write(toJSON(lt, pretty = TRUE), file = jsonFile)
  }
  stashLinkedThing(lt, type = "Txome")
}

#' @name linkedTxome
#' @rdname linkedTxome
#' 
#' @export
loadLinkedTxome <- function(jsonFile) {
  stashLinkedThing(do.call(tibble, fromJSON(jsonFile)), type="Txome")
}

### linkedTxpData -- a lightweight alternative ###

#' Make linked transcript data (linked GRanges)
#' 
#' `linkedTxpData` allows the user to save relevant _GRanges_ transcript data
#' for identifying and updating transcript metadata in a persistent manner
#' across R sessions. It can be used in combination with `inspectDigests()`
#' and `updateMetadata()`.
#' This is a lightweight version of `linkedTxome` (see
#' `makeLinkedTxome()`), which requires specifying a GTF file for building a
#' _TxDb_ and optionally a FASTA file for sequence retrieval.)
#' 
#' The `txpData` object is saved in the `getTximetaBFC()` location,
#' appended with a 32-character substring of `digest`.
#' The tibble listing all `linkedTxpData` is named
#' `linkedTxpDataTbl` and is listed in the same location.
#' 
#' @param digest character string of the full digest of the 
#' reference transcripts, see `inspectDigests()` with `fullDigest=TRUE`
#' @param indexName a name for the `index` when storing the linkedTxpData,
#' @param txpData _GRanges_ providing information about ranges 
#' representing the transcript sequences linked to `digest`
#' @param source the source of transcriptome, e.g. `denovo`. 
#' See `makeLinkedTxome()` for more information on specifying source
#' @param organism organism (e.g. "Homo sapiens")
#' @param release release number (e.g. "27")
#' @param genome genome (e.g. "GRCh38", or "none")
#' @param prefer whether to prefer loading data from linkedTxome: `txome`, 
#' linkedTxpData: `txpdata`, or the `annotated` version
#' 
#' @return nothing, the function is run for its side effects
#' 
#' @name linkedTxpData
#' @rdname linkedTxpData
#' 
#' @export
makeLinkedTxpData <- function(
  digest, 
  indexName,
  txpData,
  source,
  organism,
  release,
  genome,
  prefer = c("txome","txpdata","annotated")
) {

  message(paste0("linking user-provided metadata to digest: ",substr(digest,1,6),"..."))

  prefer <- match.arg(prefer)

  stopifnot(is(txpData, "GRanges"))

  std_sources <- c("GENCODE","Ensembl")
  source <- standardizeCapitalization(source, std_sources)
  index <- indexName

  short_digest <- substr(digest,1,6)
  digest_32 <- substr(digest,1,32)
  # a single-row tibble for the linkedTxpDataTbl
  # matches conent below in updateLinkedThingTbl()
  lt <- tibble(
    index = index,
    source = source,
    organism = organism,
    release = release,
    genome = genome,
    prefer = prefer,
    short_digest = short_digest,
    digest_32 = digest_32,
    digest = digest
  )
  
  stopifnot(nrow(lt) == 1)

  # need to save txpData in the BFC, then update the tibble
  # the name to use when saving txpData in the BFC
  # use the first 32 chars of the digest
  txpDataName <- paste0("txpdata-",digest_32)
  bfc_has <- existsInBFC(txpDataName)
  if (bfc_has) {
    message("txpData object was already saved in bfc, replacing")
    savepath <- bfcrpath(bfc, rnames=txpDataName)
  } else {
    message("saving txpData object in bfc")
    savepath <- bfcnew(bfc, txpDataName, ext=".rds")  
  }
  saveRDS(txpData, file=savepath)

  # now update the tibble to record that this exists
  stashLinkedThing(lt, type="TxpData")

}

### un-exported helper functions ###

# given a single-row tibble `lt` ("linked thing"), 
# either save this into the linkedTxomeTbl 
# or linkedTxpDataTbl, based on `type`
# (the linkedTxome/linkedTxpData tibbles lives in the tximeta BiocFileCache)
# (`type` will be used directly in message() so we use camelcase here)
stashLinkedThing <- function(lt, type=c("Txome","TxpData")) {
  type <- match.arg(type)
  tbl_name <- paste0("linked",type,"Tbl")
  stopifnot(is(lt, "tbl"))
  bfc_has_tbl <- existsInBFC(tbl_name)
  # do we need 'q' though?
  if (!bfc_has_tbl) {
    message(paste0("saving linked",type," in bfc (first time)"))
    savepath <- bfcnew(bfc, tbl_name, ext=".rds")
    linkedThingTbl <- lt
    saveRDS(linkedThingTbl, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, rnames=tbl_name)
    stopifnot(length(loadpath) == 1) # only one linkedThingTbl in BFC
    linkedThingTbl <- readRDS(loadpath)
    if (lt$index %in% linkedThingTbl$index) {
      m <- match(lt$index, linkedThingTbl$index)
      stopifnot(length(m) == 1)
      if (all(mapply(identical, lt, linkedThingTbl[m,]))) {
        message(paste0("linked",type," is same as already in bfc"))
      } else {
        message(paste0("linked",type," was different than one in bfc, replacing"))
        linkedThingTbl[m,] <- lt
      }
    } else {
      message(paste0("saving linked",type," in bfc"))
      linkedThingTbl <- rbind(linkedThingTbl, lt)
    }
    saveRDS(linkedThingTbl, file=loadpath)
  }
  invisible()
}

standardizeCapitalization <- function(source, std_sources) {
  for (src in std_sources) {
    if (tolower(source) == tolower(src)) {
      source <- src
    }
  }
  source
}

ensureColumns <- function(tbl, col_names) {
  missing_cols <- setdiff(col_names, names(tbl))
  for (m in missing_cols) {
    tbl[[m]] <- NA
  }
  tbl[col_names]
}

# this function ensures that linkedTxomes and linkedTxpData tibbles
# in the BFC have the expected columns for this version of _tximeta_
updateLinkedThingTbl <- function(type=c("Txome","TxpData")) {
  name <- paste0("linked",type,"Tbl")
  loadpath <- bfcrpath(bfc, rnames=name)
  stopifnot(length(loadpath) == 1) # only one linkedThingTbl in BFC
  linkedThingTbl <- readRDS(loadpath)
  common_cols <- c(
      "index",
      "source",
      "organism",
      "release",
      "genome",
      "prefer"
  )
  cols <- list(
    Txome = c(
      common_cols,
      "fasta",
      "gtf",
      "sha256"
    ),
    TxpData = c(
      common_cols,
      "short_digest",
      "digest_32",
      "digest"
    )
  )
  linkedThingTbl <- ensureColumns(linkedThingTbl, cols[[type]])
  saveRDS(linkedThingTbl, file=loadpath)
}

existsInBFC <- function(query) {
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, query)
  bfccount(q) > 0
}
