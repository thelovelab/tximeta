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
    loadpath <- bfcrpath(bfc, rnames="linkedTxomeTbl")
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
    loadpath <- bfcrpath(bfc, rnames=txdbName)
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
