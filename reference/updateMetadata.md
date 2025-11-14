# Update transcript metadatda for `importData()` imported data

This function takes as input a *SummarizedExperiment* as output by
[`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md),
and will update the metadata on the transcripts when possible (updating
`rowData` and/or `rowRanges` depending on the value of `ranges`).
[`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md)
uses metadata pulled from digest matches in registries used by *tximeta*
(`linkedTxome`, `linkedTxpData`, and the pre-computed digests).
Additionally, *GRanges* or *data.frame*-type data can be provided on a
one-time basis via the argument `txpData`, which will annotate
transcripts with `index="user"`. See
[`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
for how to inspect which indices have matching digests, and how to link
data to local metadata in a persistent manner.

## Usage

``` r
updateMetadata(
  se,
  txpData = NULL,
  ranges = FALSE,
  prefer = c("txome", "txpdata", "precomputed"),
  order = c("annotated", "novel", "user"),
  key = c(annotated = "tx_name", novel = "tx_name", user = "tx_name")
)
```

## Arguments

- se:

  the *SummarizedExperiment* (SE) output by
  [`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md)

- txpData:

  either *GRanges* or *data.frame*-type object to use if there is not a
  match based on digest. This is used on a one-time basis, and
  transcripts will be marked in metadata columns as
  ``` index = "user"``. See  ```makeLinkedTxome()`or`makeLinkedTxpData()\`
  for persistent metadata storage/retrieval

- ranges:

  logical, whether to add `rowRanges` (or just `rowData`)

- prefer:

  vector of length up to 3, giving the preferred order of *tximeta*'s
  transcript registries to when finding matches, with elements: `txome`:
  linkedTxome, `txpdata`: linkedTxpData, `precomputed`: the pre-computed
  digests in tximeta

- order:

  order of index, in which to update the metadata, by default the order
  is `annotation`, then `novel`, then `user`, info supplied here as
  `txpData`

- key:

  a named character vector of length 3. For each index (annotated,
  novel, and user) `key` is the name of the column to use for merging
  metadata with `rownames(se)`. The `user` index corresponds to data
  provided here as `txpData` Defaults to `"tx_name"` which often matches
  the transcript names in GENCODE

## Value

a *SummarizedExperiment* with new `rowData`, or a
*RangedSummarizedExperiment* with new metadata

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

# build custom novel GRanges data
library(GenomicRanges)
novel <- data.frame(
  seqnames = paste0("chr", rep(1:22, each=500)),  
  start = 1e6 + 1 + 0:499 * 1000, end = 1e6 + 1 + 0:499 * 1000 + 1000 - 1,
  strand = "+", tx_name = paste0("novel", 1:(22*500)),
  gene_id = paste0("novel_gene", rep(1:(22*10), each=50)), type = "protein_coding"
)
novel_gr <- as(novel, "GRanges")
names(novel_gr) <- novel$tx_name

# now update the metadata + ranges:
if (FALSE) { # \dontrun{
# this requires connection to internet (will download GENCODE GTF via FTP)
se_with_ranges <- updateMetadata(
  se, txpData=novel_gr, ranges=TRUE
)
mcols(se_with_ranges)
} # }
```
