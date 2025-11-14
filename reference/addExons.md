# Add exons to rowRanges of a transcript-level SummarizedExperiment

After running `tximeta`, the `SummarizedExperiment` output will have
`GRanges` representing the transcript locations attached as `rowRanges`
to the object. These provide the start and end of the transcript in the
genomic coordiantes, and strand information. However, the exonic
locations are not provided. This function, `addExons`, swaps out the
`GRanges` with a `GRangesList`, essentially a list along the rows of the
`SummarizedExperiment`, where each element of the list is a `GRanges`
providing the locations of the exons for that transcript.

## Usage

``` r
addExons(se)
```

## Arguments

- se:

  the SummarizedExperiment

## Value

a SummarizedExperiment

## Details

This function is designed only for transcript-level objects. This "lack
of a feature" reflects a belief on the part of the package author that
it makes more sense to think about exons belonging to transcripts than
to genes. For users desiring exonic information alongside gene-level
objects, for example, which exons are associated with a particular gene,
it is recommended to pull out the relevant `GRangesList` for the
transcripts of this gene, while the object represents transcript-level
data, such that the exons are still associated with transcripts.

For an example of `addExons`, please see the tximeta vignette.
