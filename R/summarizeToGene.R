summarizeToGene.SummarizedExperiment <- function(object,
                                                 assignRanges=c("range","abundant"),
                                                 varReduce=FALSE, ...) {

  # arguments:
  # assignRanges -
  #   "range" or "abundant", whether to return an SE
  #   with the GRanges set to the *range* of the isoforms,
  #   or the ranges of the most *abundant* isoform
  # varReduce - see ?tximport

  assignRanges <- match.arg(assignRanges)
  
  missingMetadata(object, summarize=TRUE)

  txomeInfo <- metadata(object)$txomeInfo
  txdb <- getTxDb(txomeInfo)
  message("obtaining transcript-to-gene mapping from database")

  suppressMessages({
    tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")
  })

  suppressWarnings({
    g <- getRanges(txdb=txdb, txomeInfo=txomeInfo, type="gene")
  })
      
  # need to add seqinfo for GENCODE and RefSeq
  if (all(is.na(seqlengths(g)))) {
    seqinfo(g) <- seqinfo(object)
  }
  
  txi <- list(
    abundance=assays(object)[["abundance"]],
    counts=assays(object)[["counts"]],
    length=assays(object)[["length"]],
    countsFromAbundance="no"
  )
  if ("infRep1" %in% assayNames(object)) {
    infReps <- list(assays(object)[grep("infRep", assayNames(object))])
    if (varReduce) {
      # split from per replicate list into per sample list
      # (this is what tximport expects for varReduce=TRUE)
      infReps <- list(splitInfReps(infReps[[1]]))
    }
    txi <- c(txi, infReps=infReps)
  } else {
    if (varReduce) stop("cannot calculate inferential variance without inferential replicates")
  }

  # note here about how ranges are assigned
  if (assignRanges == "abundant") {
    message("assignRanges='abundant': gene ranges assigned by isoform abundance
  see details at: ?summarizeToGene,SummarizedExperiment-method")
  } else {
    message("assignRanges='range': gene ranges assigned by total range of isoforms
  see details at: ?summarizeToGene,SummarizedExperiment-method")
  }
  
  txi.gene <- summarizeToGene(object=txi, tx2gene=tx2gene, varReduce=varReduce, ...)
  
  # put 'counts' in front to facilitate DESeqDataSet construction
  assays <- txi.gene[c("counts","abundance","length")]
  if (varReduce) {
    assays <- c(assays, txi.gene["variance"])
  } else if ("infRep1" %in% assayNames(object)) {
    infReps <- txi.gene$infReps
    assays <- c(assays, infReps)
  }

  # here do the same check/subset but with the gene-level
  # tximport assay matrices and gene ranges
  assays <- checkAssays2Txps(assays, g)

  # TODO give a warning here if there are genes in TxDb not in Salmon index?
  g <- g[rownames(assays[["counts"]])]

  # store txp IDS
  tx_ids <- CharacterList(split(tx2gene$TXNAME, tx2gene$GENEID))
  if (all(names(g) %in% names(tx_ids))) {
    tx_ids <- tx_ids[names(g)]
    mcols(g)$tx_ids <- tx_ids
  }
  
  # calculate duplication
  if ("hasDuplicate" %in% colnames(mcols(object))) {
    stopifnot(all(rownames(object) %in% tx2gene[,1]))
    t2g.o <- tx2gene[match(rownames(object),tx2gene[,1]),]
    has.dup.list <- LogicalList(split(mcols(object)$hasDuplicate, t2g.o$GENEID))
    mcols(g)$numDupObjects <- sum(has.dup.list)
  }

  # assign ranges based on most abundant isoform
  # (for the expressed genes)
  if (assignRanges == "abundant") {
    stopifnot(!is.null(rowRanges(object)))
    iso_prop <- getIsoProps(tpm=assays(object)[["abundance"]],
                            t2g=tx2gene)
    # order by the rowRanges of the outgoing object
    iso_prop <- iso_prop[names(g)]
    # remove unexpressed genes
    expressed <- !sapply(iso_prop, \(x) all(is.nan(x)))
    # compute the most abundant isoform of expressed genes
    mostAbundant <- sapply(iso_prop[expressed], \(x) names(x)[which.max(x)])
    # assign all this data to mcols(g)
    mcols(g)$isoform <- NA
    assignIdx <- names(g) %in% names(mostAbundant)
    assignNms <- names(g)[assignIdx]
    mcols(g)$isoform[assignIdx] <- mostAbundant[assignNms]
    mcols(g)$max_prop <- sapply(iso_prop, max)
    mcols(g)$iso_prop <- NumericList(iso_prop)
    # bring over new start and end positions for the subset
    # of genes that have a most abundant isoform
    txp <- rowRanges(object)
    # if addExons has been called, here we flatten
    # we could also return a GRangesList but this would be complex
    if (is(txp, "GRangesList")) {
      txp <- unlist(range(txp))
    }
    txp <- txp[mcols(g)$isoform[assignIdx]]
    # assume seqnames and strand are the same
    stopifnot(all(seqnames(g)[assignIdx] == seqnames(txp)))
    stopifnot(all(strand(g)[assignIdx] == strand(txp)))
    start(g)[assignIdx] <- start(txp)
    end(g)[assignIdx] <- end(txp)
  }
  
  metadata <- metadata(object)
  # stash countsFromAbundance value
  metadata$countsFromAbundance <- txi.gene$countsFromAbundance
  metadata$level <- "gene"
  metadata$assignRanges <- assignRanges # how were ranges assigned
  
  gse <- SummarizedExperiment(assays=assays,
                              rowRanges=g,
                              colData=colData(object),
                              metadata=metadata)  
  gse
}

#' Summarize estimated quantitites to gene-level
#'
#' Summarizes abundances, counts, lengths, (and inferential
#' replicates or variance) from transcript- to gene-level.
#' Transcript IDs are stored as a CharacterList in the \code{mcols}
#' of the output object.
#' This function operates on SummarizedExperiment objects, and
#' will automatically access the relevant TxDb (by either finding it
#' in the BiocFileCache or by building it from an ftp location).
#' This function uses the tximport package to perform summarization,
#' where a method is defined that works on simple lists.
#'
#' @param object a SummarizedExperiment produced by \code{tximeta}
#' @param assignRanges \code{"range"} or \code{"abundant"}, this argument
#' controls the way that the \code{rowRanges} of the output object are assigned
#' (note that this argument does not affect data aggregation at all).
#' The default is to just output the entire \code{range} of the gene,
#' i.e. the leftmost basepair to the rightmost basepair across all isoforms.
#' Alternatively, for expressed genes, one can obtain the start and end
#' of the most \code{abundant} isoform (averaging over all samples).
#' Non-expressed genes will have \code{range}-based positions.
#' For \code{abundant}, for expressed genes,
#' the name of the range-assigned \code{isoform}, \code{max_prop}
#' (maximum isoform proportion), and \code{iso_prop} (numeric values for
#' isoform proportions) are also returned in \code{mcols}
#' @param varReduce whether to reduce per-sample inferential replicates
#' information into a matrix of sample variances \code{variance} (default FALSE)
#' @param ... arguments passed to \code{tximport}
#'
#' @return a SummarizedExperiment with summarized quantifications
#' and transcript IDs as a CharacterList in the \code{mcols}
#'
#' @rdname summarizeToGene
#' @docType methods
#' @aliases summarizeToGene,SummarizedExperiment-method
#' 
#' @examples
#'
#' example(tximeta)
#' gse <- summarizeToGene(se)
#'
#' @importFrom tximport summarizeToGene
#' 
#' @export
setMethod("summarizeToGene", signature(object="SummarizedExperiment"),
          summarizeToGene.SummarizedExperiment)

getIsoProps <- function(tpm, t2g) {
    rownames(t2g) <- t2g[,1]
    txId <- rownames(tpm)
    stopifnot(all(txId %in% t2g[,1]))
    t2g <- t2g[txId,]
    geneId <- t2g[,2]
    tpmGene <- rowsum(tpm, geneId)
    tpmGeneExpanded <- tpmGene[match(geneId,rownames(tpmGene)),,drop=FALSE]
    isoPropMat <- tpm / tpmGeneExpanded
    isoPropMat[is.nan(isoPropMat)] <- NA
    # want to pick an isoform even if some samples have 0 expression
    split(rowMeans(isoPropMat, na.rm=TRUE), geneId)
}
