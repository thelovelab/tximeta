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
#' @param se the SummarizedExperiment
#'
#' @return a SummarizedExperiment
#'
#' @examples
#'
#' example(tximeta)
#' se <- addExons(se)
#' 
#' @export
addExons <- function(se) {

  if (metadata(se)$level == "gene") {
    stop("addExons() is design for transcript-level SummarizedExperiments, see ?addExons")
  }
  missingMetadata(se)

  txomeInfo <- metadata(se)$txomeInfo
  txdb <- getTxDb(txomeInfo)
  exons <- getRanges(txdb=txdb, txomeInfo=txomeInfo, type="exon")

  # need to add seqinfo for GENCODE and RefSeq
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
