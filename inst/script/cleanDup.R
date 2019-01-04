library(tximeta)
library(SummarizedExperiment)
files <- c("sample/quant.sf")
names <- c("foo")
coldata <- data.frame(files, names)
se <- tximeta(coldata, skipMeta=TRUE)

cleanDuplicateTxps <- function(se) {

  # ... this need to be called from within tximeta actually...
  
  # this operation requires txomeInfo
  stopifnot(!is.null(metadata(se)$txomeInfo))
  txomeInfo <- metadata(se)$txomeInfo
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
    all.dups <- dna %in% dna[dups]
    to.fix <- names(dna)[all.dups] %in% rownames(se)
    if (sum(to.fix) == 0) {
      no.dups <- TRUE
    }
  }
  if (no.dups) {
    message("no duplicated transcripts")
    return(se)
  }

  # we want to try to fix those duplicate txps in rownames(se)
  dups.to.fix <- dna[all.dups][to.fix]
  
  # get `txps`
  txdb <- getTxDb(txomeInfo)
  if (txomeInfo$source == "Ensembl") {
    suppressWarnings({
      txps <- transcripts(txdb)
    })
  } else {
    suppressWarnings({
      txps <- transcripts(txdb, columns=c("tx_id","gene_id","tx_name"))
    })      
  }
  names(txps) <- txps$tx_name

  # is there an alternative in `txps`?
  alts <- dna[!names(dna) %in% names(dups.to.fix)]
  alts <- alts[alts %in% dups.to.fix]
  alts <- alts[names(alts) %in% names(txps)]
  # no reason to prefer one since they are all in `txps`
  alts <- sort(alts[!duplicated(alts)])
  
  # only worry about dups with alternatives in `txps`
  dups.to.fix <- sort(dups.to.fix[dups.to.fix %in% alts])
  stopifnot(all(dups.to.fix == alts))

  # which rownames can we fix
  m <- match(names(dups.to.fix), rownames(se))
  # change the rownames to alternatives that are in `txps`
  rownames(se[m,]) <- names(alts)
  #...
}
