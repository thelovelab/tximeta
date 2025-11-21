# Make a DGEList from tximeta output

A simple wrapper function for constructing a DGEList for use with edgeR.
See vignette for an example. Requires installation of the edgeR package
from Bioconductor.

## Usage

``` r
makeDGEList(se, estimateDispersion = FALSE, ...)
```

## Arguments

- se:

  a SummarizedExperiment produced by tximeta

- estimateDispersion:

  logical, whether to add read-transcript ambiguity based dispersion via
  edgeR

- ...:

  arguments passed to `edgeR::estimateRTADisp`, e.g. `files` and `type`

## Value

a DGEList
