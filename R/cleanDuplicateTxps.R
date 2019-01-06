cleanDuplicateTxps <- function(missing.txps, txomeInfo, txps) {
  dup.table.id <- paste0("dups-",substr(txomeInfo$sha256,1,32))
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, dup.table.id)
  if (bfccount(q) == 1) {
    loadpath <- bfcrpath(bfc, dup.table.id)
    dup.table <- readRDS(loadpath)
  } else {
    # download the DNA from remote source
    if (is.list(txomeInfo$fasta)) {
      dna <- list()
      for (i in seq_along(txomeInfo$fasta[[1]])) {
        dna[[i]] <- readDNAStringSet(txomeInfo$fasta[[1]][i])
      }
      dna <- do.call(c, dna)
    } else {
      dna <- readDNAStringSet(txomeInfo$fasta)
    }
    # since this only concerns Ensembl txps,
    # cut off version number (and other info) from the name
    names(dna) <- sub("\\..*", "", names(dna))

    # which are duplicated based on DNA sequence
    dups <- duplicated(dna)

    no.dups <- FALSE
    if (sum(dups) == 0) {
      no.dups <- TRUE
    } else {
      # index of all seqs that participate in a dup
      all.dups <- dna %in% dna[dups]
      to.fix <- names(dna)[all.dups] %in% missing.txps
      if (sum(to.fix) == 0) {
        no.dups <- TRUE
      }
    }
    if (no.dups) {
      return(NULL)
    }

    # we want to try to fix those duplicate txps that are
    # in rownames of the assays but not in `txps`
    dups.to.fix <- dna[all.dups][to.fix]
  
    # is there an alternative in `txps`?
    alts <- dna[!names(dna) %in% names(dups.to.fix)]
    alts <- alts[alts %in% dups.to.fix]
    alts <- alts[names(alts) %in% names(txps)]
    # no reason to prefer one since they are all in `txps`
    alts <- sort(alts[!duplicated(alts)])
  
    # only worry about dups with alternatives in `txps`
    dups.to.fix <- sort(dups.to.fix[dups.to.fix %in% alts])
    stopifnot(all(dups.to.fix == alts))
    dup.table <- tibble(dups.to.fix=names(dups.to.fix),
                        alts=names(alts))

    # save it to BFC so as not to repeat the FASTA download
    savepath <- bfcnew(bfc, dup.table.id, ext=".rds")
    saveRDS(dup.table, file=savepath)
  }
  dup.table
}
