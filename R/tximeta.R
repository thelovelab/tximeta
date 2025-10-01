#' Import transcript quantification with metadata
#' 
#' The tximeta package imports abundances (TPM), estimated counts,
#' and effective lengths from quantification tools, 
#' and will output a _SummarizedExperiment_ (SE) object. 
#' For salmon and related quantification tools, [tximeta()] will
#' attempt to identify the correct provenance of the reference transcripts
#' and automatically attach the transcript ranges to the
#' SummarizedExperiment, to facilitate downstream integration with
#' other datasets. The automatic identification of reference transcripts
#' should work out-of-the-box for human or mouse transcriptomes from
#' the sources: GENCODE, Ensembl, or RefSeq. See also [tximix()] for
#' importing data when the reference transcripts were derived from 
#' a mix of annotated (e.g. GENCODE) and novel or custom transcripts.
#'
#' The main functions are:
#'   - [tximeta()] - with key argument `coldata` specifying sample information
#'   - [`summarizeToGene()`][summarizeToGene,SummarizedExperiment-method] - summarize quantification to gene-level
#'   - [tximix()] - import quantification with mixed reference transcript sets
#' 
#' All software-related questions should be posted to the Bioconductor Support Site:
#' 
#' <https://support.bioconductor.org>
#'
#' The code can be viewed at the GitHub repository,
#' which also lists the contributor code of conduct:
#'
#' <https://github.com/thelovelab/tximeta>
#' 
#' @references
#'
#'   - _tximeta_ reference:
#' 
#' Michael I. Love, Charlotte Soneson, Peter F. Hickey, Lisa K. Johnson
#' N. Tessa Pierce, Lori Shepherd, Martin Morgan, Rob Patro (2020)
#' _Tximeta: reference sequence checksums for provenance identification
#' in RNA-seq_. PLOS Computational Biology.
#' <https://doi.org/10.1371/journal.pcbi.1007664>
#'
#'   - _tximport_ reference (the effective length GLM offset and counts-from-abundance):
#' 
#' Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015)
#' _Differential analyses for RNA-seq: transcript-level estimates
#' improve gene-level inferences_. F1000Research.
#' <http://doi.org/10.12688/f1000research.7563>
#'
#' @author Michael I. Love, Charlotte Soneson, Peter Hickey, Rob Patro
#' 
#' @name tximeta-package
#' @aliases tximeta-package
#' @keywords package
"_PACKAGE"

#' Import transcript quantification with metadata
#' 
#' `tximeta` leverages the digest of the reference transcripts that were indexed
#' in order to identify metadata from the output of quantification tools. 
#' A computed digest (a hash value) can be used to uniquely identify the collection 
#' of reference sequences, and associate the dataset with other useful metadata.
#' After identification, tximeta uses a number of core Bioconductor packages (GenomicFeatures,
#' ensembldb, AnnotationHub, Seqinfo, BiocFileCache) to automatically
#' populate metadata for the user.
#' 
#' Most of the code in tximeta works to add metadata and transcript ranges
#' when the quantification was performed with salmon or related tools. However,
#' tximeta can be used with any quantification type that is supported
#' by [tximport::tximport()], where it will return an non-ranged SummarizedExperiment.
#' For other quantification tools see also the `customMetaInfo` argument below.
#' This behavior can also be triggered with `skipMeta=TRUE`.
#' 
#' tximeta performs a lookup of the digest (or hash value) of the index
#' stored in an auxilary information directory of the quantification tool's output
#' against a database of known transcriptomes, which is stored within the tximeta
#' package (`extdata/hashtable.csv`) and is continually updated to match Ensembl 
#' and GENCODE releases, with updates pushed to Bioconductor current release branch.
#' In addition, tximeta performs a lookup of the digest against a
#' locally stored table of linkedTxome references, see [makeLinkedTxome()].
#' If tximeta detects a match in either source, it will automatically populate
#' the transcript locations, the transcriptome release,
#' the genome with correct chromosome lengths, and connect the SE object to locally
#' cached derived metadata. tximeta also facilitates automatic summarization of 
#' transcript-level quantifications to the gene-level via `summarizeToGene`` without the need to 
#' manually build the correct `tx2gene` table for the reference used for indexing.
#'
#' tximeta on the first run will ask where the [BiocFileCache::BiocFileCache()] 
#' location for this package (_tximeta_) should be kept, either using a default location or a temporary
#' directory. At any point, the user can specify a location using
#' [setTximetaBFC()] and this choice will be saved for future sessions.
#' Multiple users can point to the same BiocFileCache, such that
#' transcript databases (TxDb or EnsDb) associated with certain salmon indices
#' and linkedTxomes can be accessed by different users without additional
#' effort or time spent downloading and building the relevant TxDb / EnsDb.
#' Note that, if the TxDb or EnsDb is present in AnnotationHub, tximeta will
#' use this object instead of downloading and building a TxDb/EnsDb from GTF
#' (to disable this set `useHub=FALSE`).
#'
#' In order to allow that multiple users can read and write to the
#' same location, one should set the BiocFileCache directory to
#' have group write permissions (g+w).
#'
#' @param coldata a data.frame with at least two columns (others will propogate to object):
#'   - `files` - character, paths of quantification files
#'   - `names` - character, sample names
#' if `coldata` is a vector, it is assumed to be the paths of quantification files
#' and unique sample names are created
#' @param type what quantifier was used, see [tximport::tximport()]
#' @param txOut whether to output transcript-level data.
#' `tximeta` is designed to have transcript-level output
#' with salmon, so default is `TRUE`,
#' and it's recommended to use `summarizeToGene`
#' following `tximeta` for gene-level summarization.
#' For an alevin file, `tximeta` will import the
#' gene level counts ignoring this argument (alevin
#' produces only gene-level quantification).
#' @param skipMeta whether to skip metadata generation
#' (e.g. to avoid errors if not connected to internet).
#' This calls `tximport` directly and so either
#' `txOut=TRUE` or `tx2gene` should be specified.
#' @param skipSeqinfo whether to skip the addition of Seqinfo,
#' which requires an internet connection to download the
#' relevant chromosome information table from UCSC
#' @param useHub whether to first attempt to download a TxDb/EnsDb
#' object from AnnotationHub, rather than creating from a
#' GTF file from FTP (default is TRUE). If FALSE, it will
#' force `tximeta` to download and parse the GTF
#' @param markDuplicateTxps whether to mark the status
#' (`hasDuplicate`) and names of duplicate transcripts
#' (`duplicates`) in the rowData of the SummarizedExperiment output.
#' Subsequent summarization to gene level will keep track
#' of the number of transcripts sets per gene (`numDupSets`)
#' @param cleanDuplicateTxps whether to try to clean
#' duplicate transcripts (exact sequence duplicates) by replacing
#' the transcript names that do not appear in the GTF
#' with those that do appear in the GTF
#' @param customMetaInfo the relative path to a custom metadata
#' information JSON file, relative to the paths in `files` of
#' `coldata`. For example, `customMetaInfo="meta_info.json"`
#' would indicate that in the same directory as the quantification
#' files in `files`, there are custom metadata information
#' JSON files. These should contain the SHA-256 hash of the
#' reference transcripts with the `index_seq_hash` tag
#' (see details in vignette).
#' @param skipFtp whether to avoid `ftp://` in case of
#' firewall, default is FALSE
#' @param ... arguments passed to `tximport`
#' 
#' @return a SummarizedExperiment with metadata on the `rowRanges`.
#' (if the hashed digest in the salmon or Sailfish index does not match
#' any known transcriptomes, or any locally saved `linkedTxome`,
#' `tximeta` will just return a non-ranged SummarizedExperiment)
#'
#' @examples
#'
#' # point to a salmon quantification file:
#' dir <- system.file("extdata/salmon_dm", package="tximportData")
#' files <- file.path(dir, "SRR1197474", "quant.sf") 
#' coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
#'
#' # normally we would just run the following which would download the appropriate metadata
#' # se <- tximeta(coldata)
#'
#' # for this example, we instead point to a local path where the GTF can be found
#' # by making a linkedTxome:
#' indexDir <- file.path(dir, "Dm.BDGP6.22.98_salmon-0.14.1")
#' fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
#'               "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")
#' gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
#' makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster",
#'                 release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#' se <- tximeta(coldata)
#'
#' # to clear the entire linkedTxome table
#' # (don't run unless you want to clear this table!)
#' # bfcloc <- getTximetaBFC()
#' # bfc <- BiocFileCache(bfcloc)
#' # bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays assayNames colData rowData rowRanges<- rowRanges
#' @importFrom S4Vectors metadata mcols mcols<-
#' @importFrom IRanges CharacterList LogicalList NumericList
#' @importFrom GenomicRanges seqnames strand start end start<- end<-
#' @importFrom tximport tximport summarizeToGene
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom AnnotationDbi loadDb saveDb select keys mapIds
#' @importFrom GenomicFeatures transcripts genes exonsBy cdsBy
#' @importFrom txdbmaker makeTxDbFromGFF makeTxDbFromGRanges
#' @importFrom ensembldb ensDbFromGtf EnsDb
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcnew bfcadd bfccount bfcrpath
#' @importFrom AnnotationHub AnnotationHub query dbconn dbfile
#' @importFrom Biostrings readDNAStringSet %in%
#' @importFrom tibble tibble
#' @importFrom Seqinfo Seqinfo genome<- seqlengths seqinfo seqinfo<- seqlevels
#' @importFrom tools R_user_dir file_ext
#' @importFrom utils menu packageVersion read.csv read.delim head
#' @importFrom methods is as
#'
#' @export
tximeta <- function(coldata,
                    type=NULL,
                    txOut=TRUE,
                    skipMeta=FALSE,
                    skipSeqinfo=FALSE,
                    useHub=TRUE,
                    markDuplicateTxps=FALSE,
                    cleanDuplicateTxps=FALSE,
                    customMetaInfo=NULL,
                    skipFtp=FALSE,
                    ...) {

  if (is(coldata, "vector")) {
    coldata <- data.frame(files=coldata, names=seq_along(coldata))
  }
  
  stopifnot(all(c("files","names") %in% names(coldata)))
  
  files <- as.character(coldata$files)
  names(files) <- coldata$names

  if (!all(file.exists(files))) {
    stop("the files do not exist at the location specified by 'coldata$files'")
  }

# default to salmon but print an error if files look non-salmon
  if (is.null(type)) {
    if (grepl(".quant(\\.gz)?$",coldata$files[1])) {
      stop("specify the 'type' of file to import if not salmon")
    } else {    
      type <- "salmon" # default
    }
  }

  # split out all alevin code to R/alevin.R
  # tests are in tests/testthat/test_alevin.R
  if (type == "alevin") {
    # note that `txOut` is ignored, alevin produces gene-level quantification only
    se <- tximetaAlevin(coldata = coldata, type = type, txOut = txOut,
      skipMeta = skipMeta, skipSeqinfo = skipSeqinfo, useHub = useHub, 
      markDuplicateTxps = markDuplicateTxps, cleanDuplicateTxps = cleanDuplicateTxps,
      customMetaInfo = customMetaInfo, skipFtp = skipFtp, ...)
    return(se)
  }

  message(paste("importing",type,"quantification files"))
  
  # remove the files column from colData
  coldata <- subset(coldata, select=-files)

  # metadata list with the tximeta package version, import type, and timestamp
  metadata <- makeMetadata(type)

  # when to skip attempting to load metadata
  # - skipMeta = TRUE OR 
  # - type is not a fish-method AND
  # - custom metadata file info not provided
  skipMetaLogic <- skipMeta |
    ( !type %in% c("salmon","sailfish","piscem","oarfish") &
      is.null(customMetaInfo) )
  
  if (skipMetaLogic) {
    txi <- tximport(files, type=type, txOut=txOut, ...)
    metadata$countsFromAbundance <- txi$countsFromAbundance
    se <- makeUnrangedSE(txi, coldata, metadata)
    return(se)
  } else {
    if (!txOut) stop("tximeta is designed to have transcript-level output for salmon and piscem.
  set txOut=TRUE and use summarizeToGene for gene-level summarization")
  }

  # `metaInfo` = list with quantification tool metadata from JSON files
  # either in specific directories (salmon) or alongside quantification files (newer tools)
  metaInfo <- lapply(
    files,
    getMetaInfo,
    type = type,
    customMetaInfo = customMetaInfo
  )
  
  # different styles of storing hash value by method
  hashType <- type2hashType(type)

  # check the sequence digest (hash) of the transcriptome index with 1st sample
  # readIndexSeqHash() returns a list of functions.
  # note that for oarfish, we are only looking at the `annotated_transcripts_digest`
  # for annotated + novel, use tximix...
  indexSeqHash <- readIndexSeqHash()[[hashType]](metaInfo[[1]])
  if (length(files) > 1) {
    hashes <- sapply(metaInfo, readIndexSeqHash()[[hashType]])
    if (!all(hashes == indexSeqHash)) {
      stop("the samples do not share the same index, and cannot be imported")
    }
    if (hashType == "oarfish") { 
      message("\nNote: tximeta() uses the `annotated` index digest to attach metadata,\n",
      "discarding transcripts not associated with the `annotated` index.")
      # custom check: if user is importing oarfish data and using the 'novel' flag... prompt about tximix()
      if ("novel_transcripts_digest" %in% names(metaInfo[[1]]$digest)) {
        message("\nNote: `novel` digest detected in quantification files.\n",
        "Use instead tximix(), which imports data and metadata from multiple indices.\n")
      }
    }
    checkInfReps(metaInfo)
  }

  # reshape this list object, invert the JSON hierarchy 
  # and examine consistency of the digest 'index_seq_hash'
  metaInfo <- reshapeMetaInfo(metaInfo, hashType)

  # add the per-sample metadata from quantification JSON files to the metadata list object
  metadata$quantInfo <- metaInfo
  
  # try to import files early to expose and tximport() related erreors
  txi <- tximport(files, type=type, txOut=TRUE, ...)
  metadata$countsFromAbundance <- txi$countsFromAbundance

  # use the reference seqeuence digest (hash) to try to find a match 
  # in the hash table of known and linked transcriptomes
  txomeInfo <- getTxomeInfo(indexSeqHash)
  if (is.null(txomeInfo)) {
    message("couldn't find matching transcriptome, returning non-ranged SummarizedExperiment")
    se <- makeUnrangedSE(txi, coldata, metadata)
    return(se)
  }

  # build or load a TxDb using the GTF filename as the identifier
  txdb <- getTxDb(txomeInfo, useHub=useHub, skipFtp=skipFtp)

  # build or load transcript ranges
  txps <- getRanges(txdb=txdb, txomeInfo=txomeInfo, type="txp")
  metadata$level <- "txp" # this marks the level of summarization of the SE: txp / gene

  # package up the assays from the list `txi`
  # put 'counts' in front to facilitate DESeqDataSet construction
  # and remove countsFromAbundance and infReps from assay list
  txi.nms <- c(
    "counts",
    c(setdiff(names(txi), c("counts", "countsFromAbundance", "infReps")))
  )
  assays <- txi[txi.nms]

  # if there are inferential replicates, add using rearrangeInfReps()
  if ("infReps" %in% names(txi)) {
    infReps <- rearrangeInfReps(txi$infReps)
    infReps <- addInfRepDimnames(infReps, assays)
    assays <- c(assays, infReps)
  }
  
  # Ensembl FASTA has txp version numbers,
  # but in the Ensembl GTF it is not in the txname,
  # so here we have to remove the version number to build the SummarizedExperiment
  if (txomeInfo$source %in% c("Ensembl","LocalEnsembl")) {
    txId <- sub("\\..*", "", rownames(assays[["counts"]]))
    for (nm in names(assays)) {
      rownames(assays[[nm]]) <- txId
    }
  }

  # the following function modifies assays and txps to clean duplicate txps 
  # (this occurs when salmon collapses identical transcripts during indexing)
  if (cleanDuplicateTxps) {
    dup.output.list <- duplicateTxpsClean(
      assays, txps, txomeInfo,
      markDuplicateTxps, cleanDuplicateTxps
    )
    assays <- dup.output.list$assays
    txps <- dup.output.list$txps
  }

  # special edits to rownames for GENCODE to remove chars after `|`
  # (and user didn't use --gencode when building salmon index)
  assays <- stripAllCharsAfterBar(assays)
  
  # check concordance
  assays <- checkAssays2Txps(assays, txps)
  
  # TODO we could give a warning here if there are txps in TxDb not in index
  txps <- txps[rownames(assays[["counts"]])]

  # another pass to mark duplicate transcripts
  if (markDuplicateTxps) {
    dup.output.list <- duplicateTxpsMark(
      assays, txps, txomeInfo,
      markDuplicateTxps, cleanDuplicateTxps
    )
    assays <- dup.output.list$assays
    txps <- dup.output.list$txps
  }
  
  # GENCODE and RefSeq needed Seqinfo added to seqinfo(txps)
  # function defined in `metadata_helpers.R`
  txps <- updateTxpsSeqinfo(txps, txomeInfo, skipSeqinfo)
  
  # add the txome information and TxDb information to the metadata list
  txdbInfo <- metadata(txdb)$value
  names(txdbInfo) <- metadata(txdb)$name
  metadata$txomeInfo <- txomeInfo
  metadata$txdbInfo <- txdbInfo

  se <- SummarizedExperiment(assays=assays,
                             rowRanges=txps,
                             colData=coldata,
                             metadata=metadata)
  se  
}

### un-exported functions to help tximeta() ###

# quantifiers have different location of storing index digest (hash)
type2hashType <- function(type) if (!type %in% c("piscem","oarfish")) "salmon" else type

# helper to swap across quantifiers that vary in location of the index sequence digest (hash)
readIndexSeqHash <- function() {
  list(
    salmon = function(m) m$index_seq_hash,
    piscem = function(m) m$signatures$sha256_seqs,
    oarfish = function(m) m$digest$annotated_transcripts_digest$sha256_digests$sha256_seqs
  )
}

# temporary function to map from GRCh38 to hg38 to allow easy
# comparison with UCSC objects from AnnotationHub...
# TODO we need a better solution for obtaining seqinfo for GENCODE
genome2UCSC <- function(x) {
  if (x == "GRCh38") {
    "hg38"
  } else if (x == "GRCm38") {
    "mm10"
  } else if (x == "GRCm39") {
    "mm39"
  } else {
    x
  }
}

gtf2RefSeq <- function(gtf, genome) {
  report <- sub("genomic.gff.gz","assembly_report.txt",basename(gtf))
  dir <- dirname(gtf)
  reportFtp <- paste0(dir, "/", report)
  tab <- read.delim(reportFtp, comment.char="#", header=FALSE, sep="\t", stringsAsFactors=FALSE)
  # TODO - need to figure out what to do about these un-parser friendly files
  tab <- tab[,c(7,9,10)]
  names(tab) <- c("refseqAccn","length","ucscName")
  Seqinfo(seqnames=tab$refseqAccn,
          seqlengths=tab$length,
          isCircular=NA,
          genome=genome)
}

# identify the txome based on the indexSeqHash
# - first look into the linkedTxomeTbl
# - secondly look into the pre-computed hash table in `extdata`
getTxomeInfo <- function(indexSeqHash, quiet=FALSE) {

  # first try to find any linkedTxomes in the linkedTxomeTbl
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, "linkedTxomeTbl")
  # there should only be one such entry in the tximeta bfc
  stopifnot(bfccount(q) < 2)
  if (bfccount(q) == 1) {

    # first check linkedTxomes, which should take priority over pre-computed
    loadpath <- bfcrpath(bfc, "linkedTxomeTbl")
    linkedTxomeTbl <- readRDS(loadpath)
    m <- match(indexSeqHash, linkedTxomeTbl$sha256)
    if (!is.na(m)) {
      txomeInfo <- as.list(linkedTxomeTbl[m,])
      txomeInfo$linkedTxome <- TRUE
      if (!quiet) {
        message(paste0("found matching linked transcriptome:\n[ ",
                txomeInfo$source, " - ", txomeInfo$organism,
                " - release ", txomeInfo$release," ]"))
      }
      return(txomeInfo)
      }
  }

  # if not in linkedTxomes try the pre-computed hash table...

  # TODO best this would be an external data package / future GA4GH RefGet API
  hashfile <- file.path(system.file("extdata",package="tximeta"),"hashtable.csv")
  hashtable <- read.csv(hashfile,stringsAsFactors=FALSE)
  m <- match(indexSeqHash, hashtable$sha256)
  if (!is.na(m)) {
    # now we can go get the GTF to annotate the ranges
    txomeInfo <- as.list(hashtable[m,])
    if (grepl(" ", txomeInfo$fasta)) {
      txomeInfo$fasta <- strsplit(txomeInfo$fasta, " ")
    }
    txomeInfo$linkedTxome <- FALSE
    if (!quiet) {
      message(paste0("found matching transcriptome:\n[ ",
                     txomeInfo$source, " - ", txomeInfo$organism,
                     " - release ", txomeInfo$release," ]"))
    }
    
    return(txomeInfo)
    
  }
  
  return(NULL)
}

# build or load a TxDb/EnsDb for the dataset
# useHub = whether to look in AnnotationHub for a resoruce
# skipFtp = whether to replace `ftp` with `https`
getTxDb <- function(txomeInfo, useHub=TRUE, skipFtp=FALSE) {
  # TODO what if there are multiple GTF files?
  stopifnot(length(txomeInfo$gtf) == 1)
  stopifnot(txomeInfo$gtf != "")
  txdbName <- basename(txomeInfo$gtf)
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  # look up txdbName
  q <- bfcquery(bfc, txdbName)
  # then filter for equality with rname
  q <- q[q$rname==txdbName,]

  if (skipFtp) {
    txomeInfo$gtf <- sub("ftp://","https://",txomeInfo$gtf)
  }

  ### No TxDb was found in the BiocFilecache ###
  if (bfccount(q) == 0) {

    # Ensembl and GENCODE best case we can find database on AnnotationHub
    hubSources <- c("Ensembl","GENCODE")
    srcName <- txomeInfo$source
    hubWorked <- FALSE
    if (srcName %in% hubSources) {
      ensSrc <- srcName == "Ensembl"
      dbType <- if (ensSrc) "EnsDb" else "TxDb"
      if (useHub) {
        message(paste("useHub=TRUE: checking for", dbType, "via 'AnnotationHub'"))
        ah <- AnnotationHub()
        # get records
        records <- query(ah, c(srcName, txomeInfo$organism, txomeInfo$release))
        # confirm source, organism, dbType through metadata columns
        records <- records[records$dataprovider==srcName &
                           records$species==txomeInfo$organism &
                           records$rdataclass==dbType,]        
        if (ensSrc) {
          # Confirm release number through grep on the title
          # EnsDb record titles look like "Ensembl 123 EnsDb for Homo sapiens"
          records <- records[grepl(paste(srcName, txomeInfo$release, dbType), records$title),]
        } else {
          # Narrow records based on the genome coordinates
          # GENCODE record titles look like "TxDb for Gencode v123 on hg38 coordinates"
          coords <- genome2UCSC(txomeInfo$genome)
          records <- records[grepl(coords, records$title),]
        }
        if (length(records) == 1) {
          message(paste("found matching", dbType, "via 'AnnotationHub'"))
          hubWorked <- TRUE
          txdb <- ah[[names(records)]]
          bfcadd(bfc, rname=txdbName, fpath=dbfile(dbconn(txdb)))
        } else {
          message(paste("did not find matching", dbType, "via 'AnnotationHub'"))
        }
      }
      # if check on AnnotationHub failed (or wasn't attempted)
      if (!hubWorked) {
        # build db for Ensembl
        if (ensSrc) {
          message("building EnsDb with 'ensembldb' package")
          # split code based on whether linkedTxome (bc GTF filename may be modified)
          if (!txomeInfo$linkedTxome) {
            # TODO what about suppressing all these warnings
            suppressWarnings({
              savepath <- ensDbFromGtf(
                txomeInfo$gtf,
                outfile = bfcnew(bfc, rname=txdbName, ext=".sqlite")
              )
            })
          } else {
            message("NOTE: linkedTxome with source='Ensembl', ensembldb will be used to parse GTF.
this may produce errors if the GTF is not from Ensembl, or has been modified")
            # for linkedTxome, because the GTF filename may be modified
            # we manually provide organism, genomeVersion, and version
            suppressWarnings({
              savepath <- ensDbFromGtf(
                txomeInfo$gtf,
                outfile = bfcnew(bfc, rname=txdbName, ext=".sqlite"),
                organism = txomeInfo$organism,
                genomeVersion = txomeInfo$genome,
                version = txomeInfo$release
              )
            })
          }
          txdb <- EnsDb(savepath)
        }
      }
    }

    # two cases left:
    # 1) Neither Ensembl or GENCODE source
    # 2) GENCODE source but AHub didn't work
    if ((!srcName %in% hubSources) | (srcName == "GENCODE" & !hubWorked)) {
      message("building TxDb with 'txdbmaker' package")
      # allow .rds instead of GTF
      if (tools::file_ext(txomeInfo$gtf) == "rds") {
        gtf2gr <- readRDS(txomeInfo$gtf)
        txdb <- makeTxDbFromGRanges(gtf2gr)
      } else {
        # the typical case: parse the GTF
        txdb <- makeTxDbFromGFF(txomeInfo$gtf)
      }
      saveDb(
        txdb,
        file = bfcnew(bfc, rname=txdbName, ext=".sqlite")
      )
    }

  } else {
    ### Yes, TxDb was found in the BiocFilecache ###
    loadpath <- bfcrpath(bfc, txdbName)
    if (txomeInfo$source == "Ensembl") {
      message(paste("loading existing EnsDb created:",q$create_time[1]))
      txdb <- EnsDb(loadpath)
    } else {
      message(paste("loading existing TxDb created:",q$create_time[1]))
      txdb <- loadDb(loadpath)
    }
  }
  
  txdb
}

stripAllCharsAfterBar <- function(assays) {
  testTxp <- rownames(assays[[1]])[1]
  if (grepl("ENST|ENSMUST", testTxp) & grepl("\\|", testTxp)) {
    for (i in names(assays)) {
      rownames(assays[[i]]) <- sub("\\|.*","",rownames(assays[[i]]))
    }
  }
  assays
}

# check to see if there are any missing transcripts not available
# for the rows of the tximport assay matrices. if so, give warning and subset
# (or error if all are missing)
checkAssays2Txps <- function(assays, txps) {
  assay.nms <- rownames(assays[["counts"]])
  txps.missing <- !assay.nms %in% names(txps)
  if (!all(assay.nms %in% names(txps))) {

    # it's probably ok that the messaging here uses the term 'txps',
    # because it's unlikely that we'd have genes not present in the GTF
    # which nevertheless had txps in the GTF...
    
    if (all(!assay.nms %in% names(txps))) {
      stop("none of the transcripts in the quantification files are in the GTF")
    } else {

      if (sum(txps.missing) > 3) {
        example.missing <- paste0("Example missing txps: [",
                                  paste(head(assay.nms[txps.missing],3),collapse=", "),
                                  ", ...]")
      } else {
        example.missing <- paste0("Missing txps: [",
                                  paste(assay.nms[txps.missing],collapse=", "), "]")
      }
      
      # TODO what to do here, GTF is missing some txps in FASTA for Ensembl
      warning(paste0("

Warning: the annotation is missing some transcripts that were quantified.
", sum(txps.missing), " out of ", nrow(assays[["counts"]]),
" txps were missing from GTF/GFF but were in the indexed FASTA
(e.g. this can occur with transcripts located on haplotype chromosomes).
In order to build a ranged SummarizedExperiment, these txps were removed.
To keep these txps, and to skip adding ranges, use skipMeta=TRUE

", example.missing, "
"))

      # after warning, then subset
      for (nm in names(assays)) {
        assays[[nm]] <- assays[[nm]][!txps.missing,,drop=FALSE]
      }
      
    }
  }
  assays
}

makeUnrangedSE <- function(txi, coldata, metadata) {
  assays <- txi[c("counts","abundance","length")]
  # if there are inferential replicates
  if ("infReps" %in% names(txi)) {
    infReps <- rearrangeInfReps(txi$infReps)
    assays <- c(assays, infReps)
  } else if ("variance" %in% names(txi)) {
    assays <- c(assays, txi["variance"])
  }
  if ("mean" %in% names(txi)) {
    assays <- c(assays, txi["mean"])
  }
  if ("tier" %in% names(txi)) {
    assays <- c(assays, txi["tier"])
  }
  assays <- assays[!sapply(assays, is.null)]
  SummarizedExperiment(assays=assays,
                       colData=coldata,
                       metadata=metadata)
}

# un-exported inf rep functions:

checkInfReps <- function(metaInfo) {
  if ("num_bootstraps" %in% names(metaInfo[[1]])) {
    nboot <- sapply(metaInfo, function(x) x$num_bootstraps)
    if (!all(nboot == nboot[1])) {
      message("\nNOTE: inferential replicate number not equal across files,
  may lead to errors in object construction, unless 'dropInfReps=TRUE'")
      if (any(nboot == 0)) {
        message(paste("\nNOTE: the following files (by #) have 0 inferential replicates:
  ",paste(which(nboot == 0),collapse=",")),"\n")
      }
    }
  }
}

# arrange list of inferential replicate matrices (per sample)
# into per replicate (infRep1, infRep2, ...)
rearrangeInfReps <- function(infReps) {
  nreps <- ncol(infReps[[1]])
  stopifnot(all(sapply(infReps, ncol) == nreps))
  getCols <- function(j,l) do.call(cbind, lapply(seq_along(l), function(k)  l[[k]][,j]))
  infReps <- lapply(seq_len(nreps), getCols, infReps)
  names(infReps) <- paste0("infRep",seq_len(nreps))
  infReps
}

# add dimnames from "counts" assay to list of infRep matrices
addInfRepDimnames <- function(infReps, assays) {
  lapply(infReps, function(mat) {
      rownames(mat) <- rownames(assays[["counts"]])
      colnames(mat) <- colnames(assays[["counts"]])
      mat
  })
}

# split list of inferential replicate matrices (per replicate)
# into per sample (sample1, sample2, ...)
splitInfReps <- function(infReps) {
  nsamps <- ncol(infReps[[1]])
  sample.names <- colnames(infReps[[1]])
  getCols <- function(j,l) do.call(cbind, lapply(seq_along(l), function(k)  l[[k]][,j]))
  infReps <- lapply(seq_len(nsamps), getCols, infReps)
  names(infReps) <- sample.names
  infReps
}
