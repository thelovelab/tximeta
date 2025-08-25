# code for marking or cleaning duplicate txps
duplicateTxps <- function(assays, txps, markDuplicateTxps, cleanDuplicateTxps) {
  assay.nms <- rownames(assays[["counts"]])
  txps.missing <- !assay.nms %in% names(txps)
  # either we want to mark duplicates, or clean up duplicates (if we can)
  if (markDuplicateTxps | (cleanDuplicateTxps & sum(txps.missing) > 0)) {
    dup.list <- makeDuplicateTxpsList(txomeInfo)
  }
  if (cleanDuplicateTxps & sum(txps.missing) > 0) {
    # this function swaps out rows missing in `txps`
    # for duplicate txps which are in `txps`. needed bc
    # Ensembl includes haplotype chromosome txps that duplicate
    # standard chromosome txps (identical sequence)
    missing.txps <- assay.nms[txps.missing]
    dup.table <- makeDuplicateTxpsTable(missing.txps, dup.list, names(txps))
    if (is.null(dup.table)) {
      message("no duplicated transcripts to clean")
    } else {
      message(paste("cleaning",nrow(dup.table),"duplicate transcript names"))
      # which rownames to fix
      m <- match(dup.table$dups.to.fix, assay.nms)
      stopifnot(all(!is.na(m)))
      # change the rownames to alternatives that are in `txps`
      for (nm in names(assays)) {
        assay.nms[m] <- dup.table$alts
        rownames(assays[[nm]]) <- assay.nms
      }
    }
  }

  # special edits to rownames for GENCODE to remove chars after `|`
  # (and user didn't use --gencode when building Salmon index)
  testTxp <- rownames(assays[[1]])[1]
  if (grepl("ENST|ENSMUST", testTxp) & grepl("\\|", testTxp)) {
    for (i in names(assays)) {
      rownames(assays[[i]]) <- sub("\\|.*","",rownames(assays[[i]]))
    }
  }

  assays <- checkAssays2Txps(assays, txps)
  
  # TODO we could give a warning here if there are txps in TxDb not in index
  txps <- txps[rownames(assays[["counts"]])]

  # mark duplicates in the rowData
  if (markDuplicateTxps) {
    # assay names could have changed due to cleanDuplicateTxps
    assay.nms <- rownames(assays[["counts"]])
    dups.in.rownms <- unlist(dup.list) %in% assay.nms
    dups.in.rownms <- LogicalList(split(dups.in.rownms, rep(seq_along(dup.list), lengths(dup.list))))
    names(dups.in.rownms) <- NULL
    num.dups.in.rownms <- sapply(dups.in.rownms, sum)
    just.one <- num.dups.in.rownms == 1
    if (!all(just.one)) {
      dup.list <- dup.list[just.one]
      dups.in.rownms <- dups.in.rownms[just.one]
    }
    duplicates <- dup.list[ !dups.in.rownms ]
    duplicates.id <- as.character(dup.list[ dups.in.rownms ])
    mcols(txps)$hasDuplicate <- FALSE
    mcols(txps)$duplicates <- CharacterList(as.list(rep("",length(txps))))
    if (length(duplicates) > 0) {
      message(paste(length(duplicates), "duplicate set founds"))
      mcols(txps)$hasDuplicate[ names(txps) %in% duplicates.id ] <- TRUE
      # if necessary remove any of these not in txps
      duplicates <- duplicates[ duplicates.id %in% names(txps) ]
      duplicates.id <- duplicates.id[ duplicates.id %in% names(txps) ]
      mcols(txps)$duplicates[ match(duplicates.id, names(txps)) ] <- duplicates
    } else {
      message("no duplicates found")
    }
  }
  
  return(list(assays=assays, txps=txps))
}

makeDuplicateTxpsList <- function(txomeInfo) {
  dup.list.id <- paste0("dups-",substr(txomeInfo$sha256,1,32))
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, dup.list.id)
  if (bfccount(q) == 1) {
    loadpath <- bfcrpath(bfc, dup.list.id)
    dup.list <- readRDS(loadpath)
    stopifnot(is(dup.list, "CharacterList"))
  } else {
    # download the cDNA from remote source (or load from cache)
    dna <- retrieveCDNA(txomeInfo, quiet=TRUE)
    if (txomeInfo$source == "Ensembl") {
      testTxp <- names(dna)[1]
      if (grepl("ENST|ENSMUST",testTxp)) {
        # for human and mouse, need to split on version number
        # (GTF has no version in txname, while FASTA has versions)
        names(dna) <- sub("\\..*", "", names(dna))
      } else {
        names(dna) <- sub(" .*", "", names(dna))
      }
    } else if (txomeInfo$source == "GENCODE") {
      names(dna) <- sub("\\|.*", "", names(dna))
    }
    # which are duplicated based on DNA sequence
    dups <- duplicated(dna)
    if (sum(dups) == 0) {
      dup.list <- NULL
    } else {
      # index of all seqs that participate in a dup
      all.dups <- dna %in% dna[dups]
      nms.all.dups <- names(dna)[all.dups]
      seq.all.dups <- as.character(dna[all.dups])
      dup.list <- split(nms.all.dups, seq.all.dups)
      names(dup.list) <- NULL
      dup.list <- CharacterList(dup.list)
    }
    # save it to BFC so as not to repeat the above operation
    savepath <- bfcnew(bfc, dup.list.id, ext=".rds")
    saveRDS(dup.list, file=savepath)
  }
  dup.list
}

makeDuplicateTxpsTable <- function(missing.txps, dup.list, txp.nms) {
  # dup.list is a list of the duplicate txps, according to the FASTA
  all.dups <- unlist(dup.list)

  # we want to try to fix those duplicate txps that are
  # in rownames of the assays but not in `txps`
  dups.to.fix <- intersect(all.dups, missing.txps)
  dups.to.fix.list <- LogicalList(split(all.dups %in% dups.to.fix, rep(seq_along(dup.list), lengths(dup.list))))
  dup.list <- dup.list[any(dups.to.fix.list)]
  dups.to.fix.list <- dups.to.fix.list[any(dups.to.fix.list)]
  names(dups.to.fix.list) <- NULL

  # if we have a set of duplicates where > 1 is in the rownames of assays,
  # this is unexpected (Salmon should have combined these)
  # and we pass on fixing this
  more.than.one <- sapply(dups.to.fix.list, sum) > 1
  dup.list <- dup.list[!more.than.one]
  dups.to.fix.list <- dups.to.fix.list[!more.than.one]

  # now we can name the duplicate sets by the missing txp
  names(dup.list) <- dup.list[dups.to.fix.list]
  all.dups <- unlist(dup.list)
  
  # is there an alternative in `txps`?
  dups.with.alt <- intersect(all.dups, txp.nms)
  dups.with.alt.list <- LogicalList(split(all.dups %in% dups.with.alt, rep(seq_along(dup.list), lengths(dup.list))))
  dup.list <- dup.list[any(dups.with.alt.list)]
  dups.with.alt.list <- dups.with.alt.list[any(dups.with.alt.list)]

  # pull out the alts
  names(dups.with.alt.list) <- names(dup.list)
  alts <- dup.list[dups.with.alt.list]

  # no reason to prefer one since they are all in `txps`
  alts <- sapply(alts, `[`, 1)
  stopifnot(all(alts %in% txp.nms))

  # return a table for the txp name swap
  dup.table <- data.frame(dups.to.fix=names(alts),
                          alts=unname(alts),
                          stringsAsFactors=FALSE)
  dup.table
}
