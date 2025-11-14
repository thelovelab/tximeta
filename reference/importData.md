# Import quantification across mixed reference transcripts

The *oarfish* quantification tools allows a mix of `--annotated`
reference transcripts (e.g. GENCODE, Ensembl) and `--novel` or custom
transcripts (e.g. de novo assembled transcripts not present in the
annotated set) to be used as the index for quantification.
`importData()` and associated functions facilitate import, reference
identification, and addition of metadata across `annotated` and/or
`novel` transcripts. The `importData()` function alone imports the data,
while inspection of the recognized digests and updating of transcript
metadata is handled by subsequent functions (listed in *See also*
section below).

## Usage

``` r
importData(coldata, type = "oarfish", quiet = FALSE, ...)
```

## Arguments

- coldata:

  data.frame with columns `files` and `names` as in
  [`tximeta()`](https://thelovelab.github.io/tximeta/reference/tximeta.md)

- type:

  what quantifier was used (see
  [`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html)),
  for now `importData()` works for `"oarfish"` files

- quiet:

  whether to suppress printed messages

- ...:

  arguments passed to
  [`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html)

## Value

an un-ranged SummarizedExperiment (SE) object, for use with subsequent
functions described in *See also* section

## Details

oarfish with mixed reference transcript sets may have been generated
with e.g.


    oarfish --only-index --annotated gencode.v48.transcripts.fa.gz \
      --novel my_novel_txps.fa.gz --seq-tech ont-cdna --threads 32 \
      --index-out gencode_plus_novel
    oarfish --reads reads/experiment_rep1.fastq.gz --index gencode_plus_novel \
      --output quants/experiment_rep1 --seq-tech ont-cdna \
      --filter-group no-filters --threads 32

## See also

[`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
and
[`updateMetadata()`](https://thelovelab.github.io/tximeta/reference/updateMetadata.md)
for subsequent tasks of inspecting digest matches and updating metadata,
respectively.
[`makeLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md)
can be used to add custom metadata into the registry used for inspecting
digests and then updating transcript data. A user may follow the
workflow `importData()` \>
[`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
\>
[`makeLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md)
\>
[`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
\>
[`updateMetadata()`](https://thelovelab.github.io/tximeta/reference/updateMetadata.md).
See also
[`makeLinkedTxpData()`](https://thelovelab.github.io/tximeta/reference/linkedTxpData.md)
for a lightweight alternative of linking *GRanges* metadata to a digest.

## Examples

``` r

# oarfish files using a mix of --annotated and --novel transcripts
dir <- system.file("extdata/oarfish", package="tximportData")
names <- paste0("rep", 2:4)
files <- file.path(dir, paste0("sgnex_h9_", names, ".quant.gz"))
coldata <- data.frame(files, names)

# returns an un-ranged SE object
se <- importData(coldata, type="oarfish")
#> reading in files with read.delim (install 'readr' package for speed up)
#> 1 
#> 2 
#> 3 
#> 
#> returning un-ranged SummarizedExperiment, see functions:
#> -- inspectDigests() to check matching digests
#> -- makeLinkedTxome/makeLinkedTxpData() to link digests to metadata
#> -- updateMetadata() to update metadata and optionally add ranges
```
