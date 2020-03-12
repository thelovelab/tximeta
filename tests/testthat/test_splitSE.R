context("splitSE")

test_that("splitSE works as expected", {
    ## --------------------------------------------------------------------- ##
    ## Regular dense matrix
    ## --------------------------------------------------------------------- ##
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(
            counts = matrix(1:15, nrow = 5),
            logcounts = log2(matrix(1:15, nrow = 5))
        ), 
        colData = S4Vectors::DataFrame(sID = paste0("S", 1:3),
                                       condition = c("A", "A", "B")),
        metadata = list(md1 = "annotation")
    )
    splitDf <- data.frame(spliced = c("G1", "G2", "G4"), 
                          unspliced = c("G3", "G5", "G6"),
                          stringsAsFactors = FALSE)
    ## Check that splitting fails if SE does not have row names
    expect_error(splitSE(se = se, splitDf = splitDf, assayName = "counts"),
                 regexp = "None of the IDs in the spliced, unspliced columns")
    ## Check that splitting fails if the row names of SE do not match 
    ## the IDs in splitDf
    rownames(se) <- paste0("G", c(1, 2, 33, 4, 55))
    expect_error(splitSE(se = se, splitDf = splitDf, assayName = "counts"),
                 regexp = "None of the IDs in the unspliced column")
    
    ## Set correct dimnames
    rownames(se) <- paste0("G", 1:5)
    colnames(se) <- paste0("P", 1:3)
    
    splse <- splitSE(se = se, splitDf = splitDf, assayName = "counts")
    
    expect_equal(nrow(splse), 3L)
    expect_equal(ncol(splse), 3L)
    expect_equal(colnames(splse), colnames(se))
    expect_equal(rownames(splse), splitDf[, 1])
    expect_equal(SummarizedExperiment::assayNames(splse), 
                 colnames(splitDf))
    expect_equal(S4Vectors::metadata(splse),
                 S4Vectors::metadata(se))
    expect_equal(SummarizedExperiment::colData(splse),
                 SummarizedExperiment::colData(se))
    expect_equal(assay(splse, "spliced"), 
                 matrix(c(1, 2, 4, 6, 7, 9, 11, 12, 14), 
                        nrow = 3L, ncol = 3L,
                        dimnames = list(c("G1", "G2", "G4"),
                                        c("P1", "P2", "P3"))))
    expect_equal(assay(splse, "unspliced"), 
                 matrix(c(3, 5, 0, 8, 10, 0, 13, 15, 0), 
                        nrow = 3L, ncol = 3L,
                        dimnames = list(c("G1", "G2", "G4"),
                                        c("P1", "P2", "P3"))))
    
    ## --------------------------------------------------------------------- ##
    ## Sparse matrix
    ## --------------------------------------------------------------------- ##
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(
            counts = as(matrix(1:15, nrow = 5), "sparseMatrix"),
            logcounts = log2(matrix(1:15, nrow = 5))
        ), 
        colData = S4Vectors::DataFrame(sID = paste0("S", 1:3),
                                       condition = c("A", "A", "B")),
        metadata = list(md1 = "annotation")
    )
    rownames(se) <- paste0("G", 1:5)
    colnames(se) <- paste0("P", 1:3)
    splitDf <- data.frame(spliced = c("G1", "G2", "G6"), 
                          unspliced = c("G3", "G5", "G4"),
                          stringsAsFactors = FALSE)
    
    splse <- splitSE(se = se, splitDf = splitDf, assayName = "counts")
    
    expect_equal(nrow(splse), 3L)
    expect_equal(ncol(splse), 3L)
    expect_equal(colnames(splse), colnames(se))
    expect_equal(rownames(splse), splitDf[, 1])
    expect_equal(SummarizedExperiment::assayNames(splse), 
                 colnames(splitDf))
    expect_equal(S4Vectors::metadata(splse),
                 S4Vectors::metadata(se))
    expect_equal(SummarizedExperiment::colData(splse),
                 SummarizedExperiment::colData(se))
    expect_equal(assay(splse, "spliced"), 
                 as(matrix(c(1, 2, 0, 6, 7, 0, 11, 12, 0), 
                           nrow = 3L, ncol = 3L,
                           dimnames = list(c("G1", "G2", "G6"),
                                           c("P1", "P2", "P3"))), "sparseMatrix"))
    expect_equal(assay(splse, "unspliced"), 
                 as(matrix(c(3, 5, 4, 8, 10, 9, 13, 15, 14), 
                           nrow = 3L, ncol = 3L,
                           dimnames = list(c("G1", "G2", "G6"),
                                           c("P1", "P2", "P3"))), "sparseMatrix"))
    
    ## --------------------------------------------------------------------- ##
    ## Sparse matrix, check that factors are converted to characters
    ## --------------------------------------------------------------------- ##
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(
            counts = as(matrix(1:15, nrow = 5), "sparseMatrix"),
            logcounts = log2(matrix(1:15, nrow = 5))
        ), 
        colData = S4Vectors::DataFrame(sID = paste0("S", 1:3),
                                       condition = c("A", "A", "B")),
        metadata = list(md1 = "annotation")
    )
    rownames(se) <- paste0("G", 1:5)
    colnames(se) <- paste0("P", 1:3)
    splitDf <- data.frame(spliced = c("G1", "G6", "G2"), 
                          unspliced = c("G3", "G5", "G4"),
                          stringsAsFactors = TRUE)
    
    expect_warning(splse <- splitSE(se = se, splitDf = splitDf, assayName = "counts"),
                   regexp = "Converting the 'spliced' column to character")
    
    expect_equal(nrow(splse), 3L)
    expect_equal(ncol(splse), 3L)
    expect_equal(colnames(splse), colnames(se))
    expect_equal(rownames(splse), as.character(splitDf[, 1]))
    expect_equal(SummarizedExperiment::assayNames(splse), 
                 colnames(splitDf))
    expect_equal(S4Vectors::metadata(splse),
                 S4Vectors::metadata(se))
    expect_equal(SummarizedExperiment::colData(splse),
                 SummarizedExperiment::colData(se))
    expect_equal(assay(splse, "spliced"), 
                 as(matrix(c(1, 0, 2, 6, 0, 7, 11, 0, 12), 
                           nrow = 3L, ncol = 3L,
                           dimnames = list(c("G1", "G6", "G2"),
                                           c("P1", "P2", "P3"))), "sparseMatrix"))
    expect_equal(assay(splse, "unspliced"), 
                 as(matrix(c(3, 5, 4, 8, 10, 9, 13, 15, 14), 
                           nrow = 3L, ncol = 3L,
                           dimnames = list(c("G1", "G6", "G2"),
                                           c("P1", "P2", "P3"))), "sparseMatrix"))
    
})

