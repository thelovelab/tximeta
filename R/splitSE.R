#' Split SummarizedExperiment by gene categories
#' 
#' Construct a new SummarizedExperiment by splitting one of the assays 
#' into a list of assays, each of which contains features of a given 
#' 'type'. It is assumed that there is a one-to-one correspondence 
#' between feature sets of different types; for example, these can be 
#' spliced and unspliced variants of the same transcripts. The type of
#' each feature in the original SummarizedExperiment, and the correspondence
#' between the features of different types, are given in a \code{data.frame}.
#' 
#' @param se A \code{SummarizedExperiment} object.
#' @param splitDf A \code{data.frame} with feature IDs. Each column 
#'   represents a separate feature type, and the features in a given 
#'   row are considered representatives of the same feature (and will 
#'   be represented as one feature in the output object).
#' @param assayName A character scalar, indicating the assay of \code{se}
#'   that will be split. Must be one of \code{assayNames(se)}.
#' 
#' @export
#' 
#' @return A \code{SummarizedExperiment} object with the same columns 
#'   as the input object, and the same number of assays as the number 
#'   of columns in \code{splitDf}. The assays will be named by the 
#'   column names of \code{splitDf}. The \code{colData} and \code{metadata}
#'   of the input \code{SummarizedExperiment} object are copied to the 
#'   output object. The row names are set to the feature IDs in the 
#'   first column of \code{splitDf}. 
#' 
#' @importFrom SummarizedExperiment assayNames assay SummarizedExperiment
#'   colData
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors SimpleList metadata
#'   
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = S4Vectors::SimpleList(
#'     counts = as(matrix(1:15, nrow = 5), "sparseMatrix"),
#'     logcounts = log2(matrix(1:15, nrow = 5))
#'   ), 
#'   colData = S4Vectors::DataFrame(sID = paste0("S", 1:3),
#'                                  condition = c("A", "A", "B")),
#'   metadata = list(md1 = "annotation")
#' )
#' rownames(se) <- paste0("G", 1:5)
#' colnames(se) <- paste0("P", 1:3)
#' splitDf <- data.frame(spliced = c("G1", "G2", "G6"), 
#'                       unspliced = c("G3", "G5", "G4"),
#'                       stringsAsFactors = FALSE)
#'                       
#' splse <- splitSE(se = se, splitDf = splitDf, assayName = "counts")
#' 
splitSE <- function(se, splitDf, assayName) {
    if (!(assayName %in% SummarizedExperiment::assayNames(se))) {
        stop("'assayName' must be one of assayNames(se)")
    }
    if (!is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object")
    }
    if (!is(splitDf, "data.frame")) {
        stop("'splitDf' must be a data.frame")
    }
    for (cn in colnames(splitDf)) {
        if (is(splitDf[[cn]], "factor")) {
            warning("Converting the '", cn, "' column to character")
            splitDf[[cn]] <- as.character(splitDf[[cn]])
        }
    }
    
    ## Subset splitDf to rows where at least one of the features 
    ## are in the SE.
    l <- vapply(splitDf, function(i) any(i %in% rownames(se)), FALSE)
    l <- names(l[!l])
    if (length(l) > 0) {
        stop("None of the IDs in the ", 
             paste(l, collapse = ", "), " column", ifelse(length(l) > 1, "s", ""), 
             " are present in rownames(se)")
    }
    splitDf <- splitDf[
        rowSums(sapply(splitDf, function(i) i %in% rownames(se))) > 0, 
        , drop = FALSE]
    
    ## Create new assays from the 'assayName' assay, 
    ## named by the columns of 'splitDf'.
    ## All assays will have rownames from the first 
    ## column in 'splitDf'.
    assays <- lapply(splitDf, function(i) {
        if (is(SummarizedExperiment::assay(se, assayName), "sparseMatrix")) {
            tmp <- Matrix::sparseMatrix(
                i = integer(0), j = integer(0), x = numeric(0),
                dims = c(length(i), ncol(se)),
                dimnames = list(i, colnames(se))
            )
        } else {
            tmp <- matrix(0, nrow = length(i), ncol = ncol(se),
                          dimnames = list(i, colnames(se)))
        }
        idx <- intersect(i, rownames(se))
        tmp[idx, ] <- SummarizedExperiment::assay(se, assayName)[idx, ]
        rownames(tmp) <- splitDf[, 1]
        tmp
    })
    
    SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(assays), 
        colData = SummarizedExperiment::colData(se),
        metadata = S4Vectors::metadata(se)
    )
}
