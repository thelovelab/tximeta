# Retrieve the cDNA transcript sequence for a SummarizedExperiment

This helper function retrieves the cDNA sequence of the transcripts used
for expression quantification. This function either downloads or loads
the transcript sequence from cache, it does not re-order or check
against the rows of the SummarizedExperiment (which could be already
summarized to genes for example).

## Usage

``` r
retrieveCDNA(se, quiet = FALSE)
```

## Arguments

- se:

  the SummarizedExperiment

- quiet:

  logical, suppress messages

## Value

a DNAStringSet object

## Examples

``` r

if (FALSE) { # \dontrun{
# this example is not run because it requires access to Ensembl ftp
example(tximeta)
cdna <- retrieveCDNA(se)
} # }
```
