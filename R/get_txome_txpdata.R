# identify the txome based on the digest
# - first look into the linkedTxomeTbl or linkedTxpDataTbl
# - secondly look into the pre-computed hash table in `extdata`
getTxomeInfo <- function(digest, prefer=c("txome","txpdata","precomputed"), quiet=FALSE) {
  stopifnot(all(prefer %in% c("txome","txpdata","precomputed")))
  # look through these registries in the preferred order
  for (which_registry in prefer) {
    if (which_registry %in% c("txome","txpdata")) {
      txomeInfo <- findDigestMatch(digest, which_registry, quiet)
    } else {
      txomeInfo <- findDigestMatchPrecomputed(digest, quiet)
    }
    # a non-null returned value indicates a match, return the txomeInfo
    if (!is.null(txomeInfo))
      return(txomeInfo)
  }  
  # at the end, return NULL = no match in any registry
  return(NULL)
}

findDigestMatch <- function(digest, which_registry, quiet) {
  type <- c(txome="Txome", txpdata="TxpData")[which_registry]
  stopifnot(type %in% c("Txome","TxpData"))
  bfc <- BiocFileCache(getBFCLoc())
  linkedName <- paste0("linked",type,"Tbl")
  q <- bfcquery(bfc, linkedName)
  # there should only be one such entry in the tximeta bfc
  stopifnot(bfccount(q) < 2)
  if (bfccount(q) == 1) {
    # first check linkedTxome/TxpData, which should take priority over pre-computed
    loadpath <- bfcrpath(bfc, rnames=linkedName)
    linkedTbl <- readRDS(loadpath)
    if (type == "Txome") {
      m <- match(digest, linkedTbl$sha256)
    } else if (type == "TxpData") {
      m <- match(digest, linkedTbl$digest)
    } else {
      stop("type must be one of Txome or TxpData")
    }
    if (length(m) > 1) {
      if (!quiet)
        message("found multiple matching digests, using first match")
      m <- m[1]
    }
    if (!is.na(m)) {
      txomeInfo <- as.list(linkedTbl[m,])
      txomeInfo$linkedTxome <- type == "Txome"
      txomeInfo$linkedTxpData <- type == "TxpData"
      if (!quiet) {
        message(paste0("found matching linked",type,":\n[ ",
                txomeInfo$source, " - ", txomeInfo$organism,
                " - release ", txomeInfo$release," ]"))
      }
      return(txomeInfo)
    } else {
      # no matching digest
      return(NULL)
    }
  }
  # no linked table yet
  return(NULL)
}

findDigestMatchPrecomputed <- function(digest, quiet) {
  hashfile <- file.path(system.file("extdata",package="tximeta"),"hashtable.csv")
  hashtable <- read.csv(hashfile,stringsAsFactors=FALSE)
  m <- match(digest, hashtable$sha256)
  if (!is.na(m)) {
    # now we can go get the GTF to annotate the ranges
    txomeInfo <- as.list(hashtable[m,])
    if (grepl(" ", txomeInfo$fasta)) {
      txomeInfo$fasta <- strsplit(txomeInfo$fasta, " ")
    }
    txomeInfo$linkedTxome <- FALSE
    txomeInfo$linkedTxpData <- FALSE
    if (!quiet) {
      message(paste0("found matching transcriptome:\n[ ",
                     txomeInfo$source, " - ", txomeInfo$organism,
                     " - release ", txomeInfo$release," ]"))
    }
    return(txomeInfo)
  }
  # no match in the pre-computed hash table
  return(NULL)
}

# build or load a TxDb/EnsDb for the dataset
# useHub = whether to look in AnnotationHub for a resoruce
# skipFtp = whether to replace `ftp` with `https`
getTxDb <- function(txomeInfo, useHub=TRUE, skipFtp=FALSE) {
  stopifnot(length(txomeInfo$gtf) == 1)
  stopifnot(txomeInfo$gtf != "")
  if (skipFtp) {
    txomeInfo$gtf <- sub("ftp://","https://",txomeInfo$gtf)
  }
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  txdbName <- basename(txomeInfo$gtf)
  q <- locateTxDb(txomeInfo, bfc)
  ### No TxDb was found in the BiocFilecache ###
  if (bfccount(q) == 0) {
    hubSources <- c("Ensembl","GENCODE")
    srcName <- txomeInfo$source
    hubWorked <- FALSE
    # Ensembl and GENCODE, no match in BFC
    if (srcName %in% hubSources) {
      ensSrc <- srcName == "Ensembl"
      dbType <- if (ensSrc) "EnsDb" else "TxDb"
      if (useHub) {
        txdb <- checkViaAHub(txomeInfo, srcName, ensSrc, dbType)
        hubWorked <- !is.null(txdb)
      }
      # if check on AnnotationHub failed (or wasn't attempted)
      if (!hubWorked) {
        # build db for Ensembl
        if (ensSrc) {
          txdb <- buildTxDbForEnsembl(txomeInfo, bfc, txdbName)
        }
      }
    }
    # two cases left:
    #   1. Neither Ensembl or GENCODE source
    #   2. GENCODE source but AHub didn't work
    # ...then build a TxDb with txdbmaker
    if ((!srcName %in% hubSources) | (srcName == "GENCODE" & !hubWorked)) {
      txdb <- buildTxDbFromGTF(txomeInfo, bfc, txdbName)
    }
  } else {
    ### Yes, a TxDb was found in the BiocFilecache ###
    txdb <- loadTxDbFromBFC(txomeInfo, bfc, txdbName, q)
  }
  txdb
}

locateTxDb <- function(txomeInfo, bfc) {
  txdbName <- basename(txomeInfo$gtf)
  # look up txdbName
  q <- bfcquery(bfc, txdbName)
  # filter for equality with rname
  q[q$rname==txdbName,]
}

checkViaAHub <- function(txomeInfo, srcName, ensSrc, dbType) {
  # first check for database on AnnotationHub
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
      txdb <- ah[[names(records)]]
      bfcadd(bfc, rname=txdbName, fpath=dbfile(dbconn(txdb)))
    } else {
      message(paste("did not find matching", dbType, "via 'AnnotationHub'"))
      txdb <- NULL
    }
  txdb
}

buildTxDbForEnsembl <- function(txomeInfo, bfc, txdbName) {
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
    message(paste0(
      "NOTE: linkedTxome with source='Ensembl', ensembldb will be used to parse GTF.\n",
      ". this may produce errors if the GTF is not from Ensembl, or has been modified"
    ))
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
  txdb
}

buildTxDbFromGTF <- function(txomeInfo, bfc, txdbName) {
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
  txdb
}

loadTxDbFromBFC <- function(txomeInfo, bfc, txdbName, q) {
  loadpath <- bfcrpath(bfc, rnames=txdbName)
  # check the file exists at the loadpath
  if (!file.exists(loadpath))
  stop(paste0(
    "A match was found in the digest table for this SE object,\n",
    "  but a Db not found at expected location:",
    "\n  ", loadpath
  ))
  if (txomeInfo$source == "Ensembl") {
    message(paste("loading existing EnsDb created:", q$create_time[1]))
    txdb <- ensembldb::EnsDb(loadpath)
  } else {
    message(paste("loading existing TxDb created:", q$create_time[1]))
    txdb <- AnnotationDbi::loadDb(loadpath)
  }
  txdb
}