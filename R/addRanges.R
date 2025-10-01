#' Add exons to rowRanges of a transcript-level SummarizedExperiment
#'
#' After running \code{tximeta}, the \code{SummarizedExperiment} output
#' will have \code{GRanges} representing the transcript locations
#' attached as \code{rowRanges} to the object. These provide the
#' start and end of the transcript in the genomic coordiantes, and
#' strand information. However, the exonic locations are not provided.
#' This function, \code{addExons}, swaps out the \code{GRanges}
#' with a \code{GRangesList}, essentially a list along the rows of the
#' \code{SummarizedExperiment}, where each element of the list is a
#' \code{GRanges} providing the locations of the exons for that transcript.
#'
#' This function is designed only for transcript-level objects.
#' This "lack of a feature" reflects a belief on the part of the package author
#' that it makes more sense to think about exons belonging to transcripts
#' than to genes. For users desiring exonic information alongside
#' gene-level objects, for example, which exons are associated with
#' a particular gene, it is recommended to pull out the relevant
#' \code{GRangesList} for the transcripts of this gene, while the object
#' represents transcript-level data, such that the exons are still
#' associated with transcripts.
#'
#' For an example of \code{addExons}, please see the tximeta vignette.
#' 
#' @param se the SummarizedExperiment
#'
#' @return a SummarizedExperiment
#' @export
addExons <- function(se) {

  if (metadata(se)$level == "gene") {
    stop("addExons() is design for transcript-level SummarizedExperiments, see ?addExons")
  }
  missingMetadata(se, summarize=FALSE)

  txomeInfo <- metadata(se)$txomeInfo
  txdb <- getTxDb(txomeInfo)
  exons <- getRanges(txdb=txdb, txomeInfo=txomeInfo, type="exon")

  # need to add seqinfo for (some) GENCODE and RefSeq
  if (all(is.na(seqlengths(exons)))) {
    seqinfo(exons) <- seqinfo(se)
  }
  
  # check if all transcripts are present, and then subset
  stopifnot(all(rownames(se) %in% names(exons)))
  exons <- exons[rownames(se)]
  stopifnot(all(rownames(se) == names(exons)))

  # carry over the metadata from 'se'
  mcols(exons) <- mcols(se)
  # swap the rowRanges
  rowRanges(se) <- exons
  se
}

#' Add CDS to rowRanges of a transcript-level SummarizedExperiment
#'
#' Working similarly to \code{\link{addExons}}, this function
#' can be used to add information about CDS (coding sequence)
#' to the \code{SummarizedExperiment} object. As not all transcripts
#' are coding, we have CDS information for only a subset of the
#' rows of the object. For this reason, a logical indicator for
#' whether the transcript is coding, \code{mcols(se)$coding},
#' is added as a column to the metadata columns of the \code{rowRanges}
#' of the object. An additional column, \code{mcols(se)$cds},
#' is added to the metadata columns, which is a \code{GRangesList}
#' with either the CDS regions (if the transcript is coding),
#' or the original transcript/exon ranges (if the transcript is non-coding).
#' This is necessary, as \code{GRangesList} cannot have NA elements.
#' As with \code{\link{addExons}}, this function is designed only
#' for transcript-level objects.
#' 
#' @param se the SummarizedExperiment
#'
#' @return a SummarizedExperiment
#' @export
addCDS <- function(se) {

  if (metadata(se)$level == "gene") {
    stop("addCDS() is design for transcript-level SummarizedExperiments")
  }
  missingMetadata(se, summarize=FALSE)

  txomeInfo <- metadata(se)$txomeInfo
  txdb <- getTxDb(txomeInfo)
  cds <- getRanges(txdb=txdb, txomeInfo=txomeInfo, type="cds")

  # need to add seqinfo for (some) GENCODE and RefSeq
  if (all(is.na(seqlengths(cds)))) {
    seqinfo(cds) <- seqinfo(se)
  }

  cds.in.rows <- names(cds) %in% rownames(se)
  if (!all(cds.in.rows)) {
    stopifnot(any(cds.in.rows))
    message(paste(sum(!cds.in.rows),"CDS ranges not in rows of object, dropped"))
    cds <- cds[cds.in.rows]
  }
  idx <- rownames(se) %in% names(cds)

  # note whether the transcript is coding according to CDS entry
  mcols(se)$coding <- idx
  
  nms.se <- rownames(se)[idx]
  cds <- cds[nms.se]

  # fill in with transcript or exon ranges for all transcripts
  rngs <- rowRanges(se)
  if (!is(rngs, "GRangesList")) {
    rngs <- as(rngs, "GRangesList")
  }
  mcols(se)$cds <- rngs
  stopifnot(all(rownames(se)[idx] == names(cds)))

  message("adding CDS ranges for coding transcripts, original ranges for non-coding.
see ?addCDS for more information about how CDS ranges are added")

  # add CDS ranges when the transcript is coding
  mcols(se)$cds[idx] <- cds
  
  se
}

# non-exported function:

# build or load ranges
# either transcript, exon-by-transcript, or gene ranges
getRanges <- function(txdb=txdb, txomeInfo=txomeInfo, type=c("txp","exon","cds","gene")) {
  long <- c(txp="transcript",exon="exon",cds="CDS",gene="gene")
  stopifnot(length(txomeInfo$gtf) == 1)
  stopifnot(txomeInfo$gtf != "")

  # TODO the entry in the BiocFileCache assumes that the GTF/GFF file
  # has a distinctive naming structure... works for GENCODE/Ensembl/RefSeq 
  rngsName <- paste0(type,"Rngs-",basename(txomeInfo$gtf))
  
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, rngsName)
  if (bfccount(q) == 0) {
    # now generate ranges
    message(paste("generating",long[type],"ranges"))
    # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?

    if (type == "txp") {
      ################
      ## txp ranges ##
      ################

      if (txomeInfo$source == "Ensembl") {
        suppressWarnings({
          rngs <- transcripts(txdb)
        })
      } else {
        suppressWarnings({
          rngs <- transcripts(txdb, columns=c("tx_id","gene_id","tx_name"))
        })
      }
      names(rngs) <- rngs$tx_name
      # dammit de novo transcript annotation will have
      # the transcript names as seqnames (seqid in the GFF3)
      if (tolower(txomeInfo$source) == "dammit") {
        names(rngs) <- seqnames(rngs)
      }
    } else if (type == "exon") {
      #################
      ## exon ranges ##
      #################

      # TODO suppress warnings about out-of-bound ranges for now... how to pass this on
      suppressWarnings({
        rngs <- exonsBy(txdb, by="tx", use.names=TRUE)
      })
    } else if (type == "cds") {
      #################
      ## CDS ranges ##
      #################

      # TODO suppress warnings about out-of-bound ranges for now... how to pass this on
      suppressWarnings({
        rngs <- cdsBy(txdb, by="tx", use.names=TRUE)
      })
    } else if (type == "gene") {
      #################
      ## gene ranges ##
      #################

      # TODO suppress warnings about out-of-bound ranges for now... how to pass this on
      suppressWarnings({
        rngs <- genes(txdb)
      })
    }
    savepath <- bfcnew(bfc, rngsName, ext=".rds")
    saveRDS(rngs, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, rngsName)
    message(paste("loading existing",long[type],"ranges created:",q$create_time[1]))
    rngs <- readRDS(loadpath)
  }
  rngs
}
