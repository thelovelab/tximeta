# Add CDS to rowRanges of a transcript-level SummarizedExperiment

Working similarly to
[`addExons`](https://thelovelab.github.io/tximeta/reference/addExons.md),
this function can be used to add information about CDS (coding sequence)
to the `SummarizedExperiment` object. As not all transcripts are coding,
we have CDS information for only a subset of the rows of the object. For
this reason, a logical indicator for whether the transcript is coding,
`mcols(se)$coding`, is added as a column to the metadata columns of the
`rowRanges` of the object. An additional column, `mcols(se)$cds`, is
added to the metadata columns, which is a `GRangesList` with either the
CDS regions (if the transcript is coding), or the original
transcript/exon ranges (if the transcript is non-coding). This is
necessary, as `GRangesList` cannot have NA elements. As with
[`addExons`](https://thelovelab.github.io/tximeta/reference/addExons.md),
this function is designed only for transcript-level objects.

## Usage

``` r
addCDS(se)
```

## Arguments

- se:

  the SummarizedExperiment

## Value

a SummarizedExperiment
