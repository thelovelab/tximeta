# Make a DGEList from tximeta output

A simple wrapper function for constructing a DGEList for use with edgeR.
See vignette for an example. Requires installation of the edgeR package
from Bioconductor.

## Usage

``` r
makeDGEList(se)
```

## Arguments

- se:

  a SummarizedExperiment produced by tximeta

## Value

a DGEList
