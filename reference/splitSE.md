# Split SummarizedExperiment by gene categories

Construct a new SummarizedExperiment by splitting one of the assays into
a list of assays, each of which contains features of a given 'type'. It
is assumed that there is a one-to-one correspondence between feature
sets of different types; for example, these can be spliced and unspliced
variants of the same transcripts. The type of each feature in the
original SummarizedExperiment, and the correspondence between the
features of different types, are given in a `data.frame`.

## Usage

``` r
splitSE(se, splitDf, assayName)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

- splitDf:

  A `data.frame` with feature IDs. Each column represents a separate
  feature type, and the features in a given row are considered
  representatives of the same feature (and will be represented as one
  feature in the output object).

- assayName:

  A character scalar, indicating the assay of `se` that will be split.
  Must be one of `assayNames(se)`.

## Value

A `SummarizedExperiment` object with the same columns as the input
object, and the same number of assays as the number of columns in
`splitDf`. The assays will be named by the column names of `splitDf`.
The `colData` and `metadata` of the input `SummarizedExperiment` object
are copied to the output object. The row names are set to the feature
IDs in the first column of `splitDf`.

## Examples

``` r
se <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(
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
```
