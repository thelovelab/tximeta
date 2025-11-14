# Add IDs to rowRanges of a SummarizedExperiment

For now this function just works with SummarizedExperiments with Ensembl
gene or transcript IDs. See example of usage in tximeta vignette. For
obtaining multiple matching IDs for each row of the SummarizedExperiment
set `multiVals="list"`. See `select` for documentation on use of
`multiVals`.

## Usage

``` r
addIds(se, column, fromDb = FALSE, gene = FALSE, ...)
```

## Arguments

- se:

  the SummarizedExperiment

- column:

  the name of the new ID to add (a `column` of the org package database
  or of the TxDb/EnsDb is `fromDb=TRUE`)

- fromDb:

  logical, whether to use the TxDb/EnsDb that is associated with `se`.
  Default is FALSE, and an org package is used. Currently only
  implemented for transcript level (gene=FALSE). Column names can be
  viewed with `columns(retrieveDb(se))`

- gene:

  logical, whether to map by genes or transcripts (default is FALSE). if
  rows are genes, and easily detected as such (ENSG or ENSMUSG), it will
  automatically switch to TRUE. if rows are transcripts and `gene=TRUE`,
  then it will try to use a `gene_id` column to map IDs to `column`

- ...:

  arguments passed to `mapIds`

## Value

a SummarizedExperiment

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
#> tximet> fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
#> tximet+               "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")
#> 
#> tximet> gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
#> 
#> tximet> makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster",
#> tximet+                 release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#> reading digest from indexDir: .../Dm.BDGP6.22.98_salmon-0.14.1
#> NOTE: this digest matches one in the pre-computed digest table
#> saving linkedTxome in bfc (first time)
#> 
#> tximet> se <- tximeta(coldata)
#> importing salmon quantification files
#> reading in files with read.delim (install 'readr' package for speed up)
#> 1 
#> 
#> found matching linkedTxome:
#> [ LocalEnsembl - Drosophila melanogaster - release 98 ]
#> building TxDb with 'txdbmaker' package
#> Import genomic features from the file as a GRanges object ... 
#> OK
#> Prepare the 'metadata' data frame ... 
#> OK
#> Make the TxDb object ... 
#> Warning: genome version information is not available for this TxDb object
#> OK
#> generating transcript ranges
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
library(org.Dm.eg.db)  
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     I, expand.grid, unname
#> 
se <- addIds(se, "REFSEQ", gene=FALSE)
#> mapping to new IDs using org.Dm.eg.db
#> if all matching IDs are desired, and '1:many mappings' are reported,
#> set multiVals='list' to obtain all the matching IDs
#> 'select()' returned 1:many mapping between keys and columns
```
