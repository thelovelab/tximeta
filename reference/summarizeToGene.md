# Summarize estimated quantitites to gene-level

Summarizes abundances, counts, lengths, (and inferential replicates or
variance) from transcript- to gene-level. Transcript IDs are stored as a
CharacterList in the `mcols` of the output object. This function
operates on SummarizedExperiment objects, and will automatically access
the relevant TxDb (by either finding it in the BiocFileCache or by
building it from an ftp location). This function uses the tximport
package to perform summarization, where a method is defined that works
on simple lists.

## Usage

``` r
# S4 method for class 'SummarizedExperiment'
summarizeToGene(
  object,
  assignRanges = c("range", "abundant"),
  varReduce = FALSE,
  skipRanges = FALSE,
  ...
)
```

## Arguments

- object:

  a SummarizedExperiment produced by `tximeta`

- assignRanges:

  `"range"` or `"abundant"`, this argument controls the way that the
  `rowRanges` of the output object are assigned (note that this argument
  does not affect data aggregation at all). The default is to just
  output the entire `range` of the gene, i.e. the leftmost basepair to
  the rightmost basepair across all isoforms. Alternatively, for
  expressed genes, one can obtain the start and end of the most
  `abundant` isoform (averaging over all samples). Non-expressed genes
  will have `range`-based positions. For `abundant`, for expressed
  genes, the name of the range-assigned `isoform`, `max_prop` (maximum
  isoform proportion), and `iso_prop` (numeric values for isoform
  proportions) are also returned in `mcols`

- varReduce:

  whether to reduce per-sample inferential replicates information into a
  matrix of sample variances `variance` (default FALSE)

- skipRanges:

  whether to skip making use of, or outputting, rowRanges (default
  FALSE)

- ...:

  arguments passed to `tximport`

## Value

a SummarizedExperiment with summarized quantifications and transcript
IDs as a CharacterList in the `mcols`

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
#> loading existing TxDb created: 2026-04-13 17:26:46
#> loading existing transcript ranges created: 2026-04-13 17:26:47
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
gse <- summarizeToGene(se)
#> loading existing TxDb created: 2026-04-13 17:26:46
#> obtaining transcript-to-gene mapping from database
#> generating gene ranges
#> assignRanges='range': gene ranges assigned by total range of isoforms
#>   see details at: ?summarizeToGene,SummarizedExperiment-method
#> summarizing abundance
#> summarizing counts
#> summarizing length
```
