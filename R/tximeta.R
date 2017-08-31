#' tximeta
#'
#' tximeta is magic
#'
#' @param coldata a data.frame with at least two columns:
#' \itemize{
#' \item{"files"}{character vector pointing to sample quantification files}
#' \item{"names"}{character vector of sample names}
#' }
#' 
#' @return a SummarizedExperiment
#'
#' @import SummarizedExperiment
#' @importFrom tximport tximport
#' @importFrom rjson fromJSON
#' @importFrom GenomicFeatures makeTxDbFromGFF
#'
#' @export
tximeta <- function(coldata) {
  stopifnot(all(c("files","names") %in% names(coldata)))
  files <- coldata$files
  names(files) <- coldata$names
  dir <- dirname(files)
  jsonPath <- file.path(dir, "aux_info/meta_info.json")
  stopifnot(file.exists(jsonPath))
  metaInfo <- fromJSON(file=jsonPath[1])
  indexSeqHash <- metaInfo$index_seq_hash
  hashtable <- read.csv(here("extdata","hashtable.csv"),stringsAsFactors=FALSE)
  idx <- match(indexSeqHash, hashtable$index_seq_hash)
  stopifnot(length(idx) == 1)
  message(with(hashtable[idx,,drop=FALSE],
               paste0("found matching transcriptome:\n[ ",
                      source," - ",organism," - version ",version," ]")))
  gtf <- hashtable$gtf[idx]
  txdbName <- paste0(basename(gtf),".sqlite")
  if (!file.exists(txdbName)) {
    message("building transcript database")
    txdb <- makeTxDbFromGFF(gtf)
    saveDb(txdb, file=txdbName)
  } else {
    message("loading existing transcript database")
    txdb <- loadDb(txdbName)
  }
  message("generating transcript GRanges")
  txps <- transcripts(txdb)

  # TODO fix seqinfo(txps)
  
  names(txps) <- txps$tx_name
  txi <- tximport(files, type="salmon", txOut=TRUE)
  stopifnot(all(rownames(txi$abundance) %in% names(txps)))
  txps <- txps[rownames(txi$abundance)]
  coldataSub <- subset(coldata, select=-files)

  # prepare metadata list
  tximetaInfo <- list(version=packageVersion("tximeta"),
                      importTime=Sys.time())
  
  metadata <- list(
    quantInfo=metaInfo,
    txomeInfo=as.list(hashtable[idx,]),
    tximetaInfo=tximetaInfo
  )
  
  se <- SummarizedExperiment(assays=txi[c("abundance","counts","length")],
                             rowRanges=txps,
                             colData=coldataSub,
                             metadata=metadata)
  se
  
}
