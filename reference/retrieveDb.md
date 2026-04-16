# Retrieve the TxDb or EnsDb associated with a SummarizedExperiment

SummarizedExperiment objects returned by
[`tximeta`](https://thelovelab.github.io/tximeta/reference/tximeta.md)
have associated TxDb or EnsDb databases which are cached locally and
used to perform various metadata related tasks. This helper function
retrieves the database itself for the user to perform any additional
operations.

## Usage

``` r
retrieveDb(se)
```

## Arguments

- se:

  the SummarizedExperiment

## Value

a database object

## Examples

``` r

example(tximeta)
#> 
#> tximet> # point to a salmon quantification file:
#> tximet> dir <- system.file("extdata/salmon_dm", package="tximportData")
#> 
#> tximet> files <- file.path(dir, "SRR1197474", "quant.sf") 
#> 
#> tximet> coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
#> 
#> tximet> # normally we would just run the following which would download the appropriate metadata
#> tximet> # se <- tximeta(coldata)
#> tximet> 
#> tximet> # for this example, we instead point to a local path where the GTF can be found
#> tximet> # by making a linkedTxome:
#> tximet> indexDir <- file.path(dir, "Dm.BDGP6.22.98_salmon-0.14.1")
#> 
#> tximet> dmFTP <- "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/"
#> 
#> tximet> fastaFTP <- paste0(
#> tximet+   dmFTP,
#> tximet+   c("cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
#> tximet+     "ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")
#> tximet+ )
#> 
#> tximet> gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
#> 
#> tximet> makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster",
#> tximet+                 release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#> reading digest from indexDir: .../Dm.BDGP6.22.98_salmon-0.14.1
#> NOTE: this digest matches one in the pre-computed digest table
#> linkedTxome metadata was same as already in bfc
#> 
#> tximet> se <- tximeta(coldata)
#> importing salmon quantification files
#> reading in files with read.delim (install 'readr' package for speed up)
#> 1 
#> 
#> found matching linkedTxome:
#> [ LocalEnsembl - Drosophila melanogaster - release 98 ]
#> loading existing TxDb created: 2026-04-16 14:26:19
#> Loading required package: GenomicFeatures
#> Loading required package: Seqinfo
#> Loading required package: GenomicRanges
#> loading existing transcript ranges created: 2026-04-16 14:26:20
#> Warning: 
#> 
#> Warning: the annotation is missing some transcripts that were quantified.
#> 5 out of 33706 txps were missing from GTF/GFF but were in the indexed FASTA
#> (e.g. this can occur with transcripts located on haplotype chromosomes).
#> In order to build a ranged SummarizedExperiment, these txps were removed.
#> To keep these txps, and to skip adding ranges, use skipMeta=TRUE
#> 
#> Example missing txps: [FBtr0307759, FBtr0084079, FBtr0084080, ...]
#> 
#> tximet> # to clear the entire linkedTxome table
#> tximet> # (don't run unless you want to clear this table!)
#> tximet> # bfcloc <- getTximetaBFC()
#> tximet> # bfc <- BiocFileCache(bfcloc)
#> tximet> # bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
#> tximet> 
#> tximet> 
#> tximet> 
#> tximet> 
edb <- retrieveDb(se)
#> loading existing TxDb created: 2026-04-16 14:26:19
```
