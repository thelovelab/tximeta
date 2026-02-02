# this initializes the metadata list with basic info
makeMetadata <- function(type) {
  tximetaInfo <- list(version=packageVersion("tximeta"),
                      type=type,
                      importTime=Sys.time())
  list(tximetaInfo=tximetaInfo)
}

# read metadata files from Salmon/piscem/oarfish output
# customMetaInfo = path of the custom metadata info file
getMetaInfo <- function(file, type, customMetaInfo=NULL) {
  dir <- dirname(file)

  # users can specify any arbitrary location for the metadata,
  # allowing for any quantification tool to be paired with tximeta.
  # we first deal with this case, then move to Salmon and piscem
  if (!is.null(customMetaInfo)) {
    jsonPath <- file.path(dir, customMetaInfo)

    # salmon or piscem have different metadata locations,
    # so we handle these separately...
  } else {

    # salmon:
    if (type == "salmon") {
      # the default Salmon auxiliary information location
      auxDir <- "aux_info" 
      if (!file.exists(file.path(dir, auxDir))) {
        auxDir <- customAuxDir(dir, auxDir)
      }
      # read in the metadata
      jsonPath <- file.path(dir, auxDir, "meta_info.json")

      # piscem / oarfish:
    } else if (type %in% c("piscem","oarfish")) {

      # read in the metadata
      quantFile <- basename(file)
      metadataFile <- sub(".quant(\\.gz)?$", ".meta_info.json", quantFile)
      jsonPath <- file.path(dir, metadataFile)
      
    } else {
      stop("expected type = 'salmon', 'piscem', 'oarfish'")
    }
  }
  if (!file.exists(jsonPath)) {
    stop("\n\n  the quantification files exist, but the metadata files are missing.
  tximeta (and other downstream software) require the entire output directory
  of Salmon/alevin, or for piscem/oarfish the metadata files to be colocated with the
  quant files. The total output of Salmon/alevin/piscem/oarfish includes files with
  critical metadata for tximeta to work. Alternatively, you can set
  skipMeta=TRUE or use tximport \n\n") 
  }
  fromJSON(jsonPath)
}

# Salmon allows users to change the name of the auxiliary directory
# just in case this was changed by the user...
customAuxDir <- function(dir, auxDir) {
  jsonPath <- file.path(dir, "cmd_info.json")
  if (!file.exists(jsonPath)) {
    stop("metadata files are missing, tximeta requires the full Salmon/piscem output files")
  }
  cmd_info <- jsonlite::fromJSON(jsonPath)
  if ("auxDir" %in% names(cmd_info)) {
    auxDir <- cmd_info$auxDir
  }
  auxDir
}

# reshape metadata info from salmon or other quantifiers
reshapeMetaInfo <- function(metaInfo, hashType) {
  unionTags <- unique(unlist(lapply(metaInfo, names)))
  # re-order by tag t, then sample i
  out <- lapply(unionTags, function(t) {
    sapply(seq_along(metaInfo), function(i) {
      metaInfo[[i]][[t]]
    })
  })
  names(out) <- unionTags
  if (all(out$eq_class_properties == list())) {
    out$eq_class_properties <- NULL
  }
  # TODO have alternate testing for piscem or oarfish?
  if (hashType == "salmon") {
    stopifnot(all(out$index_seq_hash == out$index_seq_hash[1]))
    stopifnot(all(out$index_name_hash == out$index_name_hash[1]))
    out$index_seq_hash <- out$index_seq_hash[1]
    out$index_name_hash <- out$index_name_hash[1]
  }
  out
}

updateTxpsSeqinfo <- function(txps, txomeInfo, skipSeqinfo) {
  # Ensembl already has nice seqinfo attached, nothing needed

  missingSeqinfo <- any(is.na(seqlengths(txps)))

  # if GENCODE, and not from AHub (which have seqinfo)...
  if (txomeInfo$source == "GENCODE" & !skipSeqinfo & missingSeqinfo) {
    message("fetching genome info for GENCODE")
    ucsc.genome <- genome2UCSC(txomeInfo$genome)
    try(seqinfo(txps) <- Seqinfo(genome=ucsc.genome)[seqlevels(txps)])
  } else if (txomeInfo$source == "RefSeq" & !skipSeqinfo & missingSeqinfo) {
    
    # if RefSeq...
    message("fetching genome info for RefSeq")
    refseq.genome <- gtf2RefSeq(txomeInfo$gtf, txomeInfo$genome)
    stopifnot(all(seqlevels(txps) %in% seqnames(refseq.genome)))
    try(seqinfo(txps) <- refseq.genome[seqlevels(txps)])
  }

  txps
}

missingMetadata <- function(se, summarize=FALSE) {
  msg <- "transcriptome metadata is missing: `metadata(se)$txomeInfo`
  either: 
    
  (1) the object was not produced by tximeta.
      Note that importData() is not supported yet, 
      see updateMetadata() instead, or
  (2) tximeta could not recognize the digest of the transcriptome.
  
  If (2), use a linkedTxome to provide the missing metadata and re-run tximeta()"

  if (summarize) {
    msg <- paste0(msg, "
  or provide a `tx2gene` data.frame and set `skipRanges=TRUE`")
  }

  if (is.null(metadata(se)$txomeInfo)) stop(msg)
  
}
