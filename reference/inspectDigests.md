# Inspect digest matches from `importData()` imported data

This function takes as input a *SummarizedExperiment* as output by
[`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md)
and returns a tibble with information about the digest-match status of
two indices (`annotated` and `novel`), with respect to *tximeta*
metadata. Inspection of index digests can be run iteratively, checking
if the digests used in the mixed reference transcript set have a match
against 1) pre-computed digests representing standard annotated sets
(e.g. GENCODE, Ensembl, etc.) or 2) digests added by the user to a local
registry with
[`makeLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md)
(GTF file) or `makeLinkedTxpData` (*GRanges*-based metadata). Optional
columns may be added if specified by `fullDigest=TRUE` (include the full
digest) and/or `count=TRUE` (add matching transcript ID counts per
index). Following inspection, one can run
[`updateMetadata()`](https://thelovelab.github.io/tximeta/reference/updateMetadata.md)
to automatically update the transcript metadata using the sources
indicated by this function.

## Usage

``` r
inspectDigests(
  se,
  type = "oarfish",
  prefer = c("txome", "txpdata", "precomputed"),
  fullDigest = FALSE,
  count = FALSE
)
```

## Arguments

- se:

  the *SummarizedExperiment* output by
  [`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md),
  or alternatively just `metadata(se)$quantInfo`, a list of metadata
  information from the quantification tool (assuming `annotated` and
  `novel` indices both used)

- type:

  what quantifier was used (see
  [`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html))

- prefer:

  vector of length up to 3, giving the preferred order of *tximeta*'s
  transcript registries to when finding matches, with elements: `txome`:
  linkedTxome, `txpdata`: linkedTxpData, `precomputed`: the pre-computed
  digests in tximeta

- fullDigest:

  whether to include the full digest string in the output, in addition
  to the shortened 6-char version

- count:

  whether to count the number of matching transcripts ID to each index
  (only possible for those indices that have matching metadata).
  Counting requires loading transcript data, either from locally cached
  databases or from GTF files.

## Value

a 2-row tibble of the `annotated` and `novel` index, their matching
information if available (source, organism, release), for matches,
whether it is a `linkedTxome` or a `linkedTxpData` (both \`FALSE“ for
pre-computed) and a small 6 character version of the digest itself.

## Examples

``` r

example(importData)
#> 
#> imprtD> # oarfish files using a mix of --annotated and --novel transcripts
#> imprtD> dir <- system.file("extdata/oarfish", package="tximportData")
#> 
#> imprtD> names <- paste0("rep", 2:4)
#> 
#> imprtD> files <- file.path(dir, paste0("sgnex_h9_", names, ".quant.gz"))
#> 
#> imprtD> coldata <- data.frame(files, names)
#> 
#> imprtD> # returns an un-ranged SE object
#> imprtD> se <- importData(coldata, type="oarfish")
#> reading in files with read.delim (install 'readr' package for speed up)
#> 1 
#> 2 
#> 3 
#> 
#> returning un-ranged SummarizedExperiment, see functions:
#> -- inspectDigests() to check matching digests
#> -- makeLinkedTxome/makeLinkedTxpData() to link digests to metadata
#> -- updateMetadata() to update metadata and optionally add ranges
# now we have an `se` created by importData()...
inspectDigests(se)
#> # A tibble: 2 × 8
#>   index    source organism release genome linkedTxome linkedTxpData small_digest
#>   <chr>    <chr>  <chr>    <chr>   <chr>  <lgl>       <lgl>         <chr>       
#> 1 annotat… GENCO… Homo sa… 48      GRCh38 FALSE       FALSE         6fc626      
#> 2 novel    NA     NA       NA      NA     NA          NA            43158f      
# can then update the registry via makeLinkedTxome() and re-run inspection
```
