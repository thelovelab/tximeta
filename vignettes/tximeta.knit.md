---
title: "tximeta: transcript quantification import with automatic metadata"
author: "Michael Love, Rob Patro, Charlotte Soneson, Peter Hickey"
date: "07/16/2018"
output: 
  rmarkdown::html_document:
    highlight: tango
abstract: >
  `tximeta` performs numerous annotation and metadata gathering tasks on
  behalf of users during the import of transcript quantifications from
  *Salmon* or *Sailfish* into R/Bioconductor. Transcript ranges and
  metadata are added automatically, facilitating combining multiple
  genomic datasets and helping to prevent bioinformatic errors.
vignette: |
  %\VignetteIndexEntry{Transcript quantification import with automatic metadata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Analysis starts with sample table

The first step using `tximeta` is to read in the sample table, which
will become the *column data*, `colData`, of the final object,
a *SummarizedExperiment*. The sample table should contain all the
information we need to identify the *Salmon* quantification
directories. Here we will use a *Salmon* quantification file in the
*tximportData* package to demonstrate the usage of `tximeta`. We do
not have a sample table, so we construct one in R. It is recommended
to keep a sample table as a CSV or TSV file while working on an
RNA-seq project with multiple samples.


```r
dir <- system.file("extdata/salmon_dm", package="tximportData")
# here gzipped, normally these are not
files <- file.path(dir, "SRR1197474_cdna", "quant.sf.gz") 
file.exists(files)
```

```
## [1] TRUE
```

```r
coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
coldata
```

```
##                                                                                          files
## 1 /Users/love/Library/R/3.5/library/tximportData/extdata/salmon_dm/SRR1197474_cdna/quant.sf.gz
##        names condition
## 1 SRR1197474         A
```

`tximeta` expects at least two columns in `coldata`: 

1. `files` - a pointer to the `quant.sf` files
2. `names` - the unique names that should be used to identify samples

# Running tximeta from a sample table

Normally, we would just run `tximeta` like so:


```r
library(tximeta)
se <- tximeta(coldata)
```

However, to avoid downloading remote GTF files during this vignette, we
will point to a GTF file saved locally (in the *tximportData*
package). We link the transcriptome of the *Salmon* index to its
locally saved GTF. 

This following code is therefore not recommended
for a typically workflow, but is particular to the vignette code.


```r
dir <- system.file("extdata", package="tximeta")
indexDir <- file.path(dir, "Drosophila_melanogaster.BDGP6.cdna.v92_salmon_0.10.2")
fastaFTP <- "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
dir2 <- system.file("extdata/salmon_dm", package="tximportData")
gtfPath <- file.path(dir2,"Drosophila_melanogaster.BDGP6.92.gtf.gz")
makeLinkedTxome(indexDir=indexDir,
                source="Ensembl",
                organism="Drosophila melanogaster",
                version="92",
                genome="BDGP6",
                fasta=fastaFTP,
                gtf=gtfPath,
                write=FALSE)
```

```
## saving linkedTxome in bfc (first time)
```


```r
library(tximeta)
se <- tximeta(coldata)
```

```
## importing quantifications
```

```
## reading in files with read_tsv
```

```
## 1 
## found matching linked transcriptome:
## [ Ensembl - Drosophila melanogaster - version 92 ]
## loading existing EnsDb created: 2018-07-14 14:45:15
## generating transcript ranges
```

```
## Warning in checkTxi2Txps(txi, txps): missing some transcripts!
##  9 out of 30061 are missing from the GTF and dropped from SummarizedExperiment output
```

# What happened? 

`tximeta` recognized the signature of the transcriptome that the files
were quantified against, it accessed the GTF file of the transcriptome
source, found and attached the transcript ranges, and added the
appropriate transcriptome and genome metadata.  A remote GTF is only
downloaded once, and a local or remote GTF is only parsed to build a
*TxDb* once: if `tximeta` recognizes that it has seen this *Salmon*
index before, it will use a cached version of the transcript metadata
and ranges.

Note the warning above that 9 of the transcripts are missing from the
GTF file and so are dropped from the final output. This is a problem
coming from the annotation source, and not easily avoided by
`tximeta`. 

We plan to create and maintain a large table of signatures for as many
sources, organisms, versions of transcriptomes as possible. We are
also developing support for *linked transcriptomes*, where one or
more sources for transcript sequences have been combined or
filtered. See the `linkedTxome` vignette in this package for a
demonstration. (The *makeLinkedTxome* function was used above to avoid
downloading the GTF during the vignette building process.)

# Examining SummarizedExperiment output

We, of course, have our coldata from before. Note that we've removed `files`.


```r
colData(se)
```

```
## DataFrame with 1 row and 2 columns
##                  names   condition
##            <character> <character>
## SRR1197474  SRR1197474           A
```

Here we show the three matrices that were imported. (Note: this part
would need updating for un-reduced inferential variance matrices.) 
(Second note: the downstream packages would need updating of their
functions or workflows, e.g. `DESeqDataSetFromTximport` needs a little
update to work with these *SummarizedExperiments* instead of simple lists.)


```r
assayNames(se)
```

```
## [1] "counts"    "abundance" "length"
```

`tximeta` has imported the correct ranges for the transcripts:


```r
rowRanges(se)
```

```
## GRanges object with 30052 ranges and 6 metadata columns:
##               seqnames            ranges strand |       tx_id
##                  <Rle>         <IRanges>  <Rle> | <character>
##   FBtr0075502       3L 15808322-15808883      + | FBtr0075502
##   FBtr0300738       3R   5783105-5787336      + | FBtr0300738
##   FBtr0300739       3R   5781762-5787336      + | FBtr0300739
##   FBtr0300737       3R   5781762-5787336      + | FBtr0300737
##   FBtr0300736       3R   5783105-5787336      + | FBtr0300736
##           ...      ...               ...    ... .         ...
##   FBtr0347432        Y   3074924-3075180      + | FBtr0347432
##   FBtr0307579        X 21156259-21156621      - | FBtr0307579
##   FBtr0089614       3R 30212903-30213142      + | FBtr0089614
##   FBtr0299927        X   7897250-7897987      - | FBtr0299927
##   FBtr0303313       3R 12893972-12894529      - | FBtr0303313
##                   tx_biotype tx_cds_seq_start tx_cds_seq_end     gene_id
##                  <character>        <integer>      <integer> <character>
##   FBtr0075502 protein_coding         15808418       15808716 FBgn0036531
##   FBtr0300738 protein_coding          5783217        5787117 FBgn0037375
##   FBtr0300739 protein_coding          5781900        5787117 FBgn0037375
##   FBtr0300737 protein_coding          5781900        5787117 FBgn0037375
##   FBtr0300736 protein_coding          5783217        5787117 FBgn0037375
##           ...            ...              ...            ...         ...
##   FBtr0347432     pseudogene             <NA>           <NA> FBgn0267873
##   FBtr0307579     pseudogene             <NA>           <NA> FBgn0052511
##   FBtr0089614     pseudogene             <NA>           <NA> FBgn0000281
##   FBtr0299927     pseudogene             <NA>           <NA> FBgn0260447
##   FBtr0303313     pseudogene             <NA>           <NA> FBgn0053929
##                   tx_name
##               <character>
##   FBtr0075502 FBtr0075502
##   FBtr0300738 FBtr0300738
##   FBtr0300739 FBtr0300739
##   FBtr0300737 FBtr0300737
##   FBtr0300736 FBtr0300736
##           ...         ...
##   FBtr0347432 FBtr0347432
##   FBtr0307579 FBtr0307579
##   FBtr0089614 FBtr0089614
##   FBtr0299927 FBtr0299927
##   FBtr0303313 FBtr0303313
##   -------
##   seqinfo: 25 sequences from BDGP6 genome
```

We have appropriate genome information, which prevents us from making 
bioinformatic mistakes:


```r
seqinfo(se)
```

```
## Seqinfo object with 25 sequences from BDGP6 genome:
##   seqnames             seqlengths isCircular genome
##   211000022278279           12714       <NA>  BDGP6
##   211000022278436            2815       <NA>  BDGP6
##   211000022278449            1947       <NA>  BDGP6
##   211000022278760            1144       <NA>  BDGP6
##   211000022279165            1118       <NA>  BDGP6
##   ...                         ...        ...    ...
##   Unmapped_Scaffold_8       88768       <NA>  BDGP6
##   X                      23542271       <NA>  BDGP6
##   Y                       3667352       <NA>  BDGP6
##   mitochondrion_genome      19517       <NA>  BDGP6
##   rDNA                      76973       <NA>  BDGP6
```

# Easy summarization to gene-level

Because the SummarizedExperiment maintains all the metadata of its
creation, it also keeps a pointer to the necessary database for
summarizing transcript-level quantifications and bias corrections to
the gene-level. If necessary, `summarizeToGene` can pull down the
remote source for summarization, but given that we've already built a
TxDb once, it simply loads the stashed version.


```r
gse <- summarizeToGene(se)
```

```
## loading existing EnsDb created: 2018-07-14 14:45:15
```

```
## obtaining transcript-to-gene mapping from TxDb
```

```
## summarizing abundance
```

```
## summarizing counts
```

```
## summarizing length
```

```r
rowRanges(gse)
```

```
## GRanges object with 14026 ranges and 6 metadata columns:
##               seqnames            ranges strand |     gene_id   gene_name
##                  <Rle>         <IRanges>  <Rle> | <character> <character>
##   FBgn0000008       2R 22136968-22172834      + | FBgn0000008           a
##   FBgn0000014       3R 16807214-16830049      - | FBgn0000014       abd-A
##   FBgn0000015       3R 16927212-16972236      - | FBgn0000015       Abd-B
##   FBgn0000017       3L 16615866-16647882      - | FBgn0000017         Abl
##   FBgn0000018       2L 10973443-10975293      - | FBgn0000018         abo
##           ...      ...               ...    ... .         ...         ...
##   FBgn0285958       3L 11518798-11522713      - | FBgn0285958        Fuca
##   FBgn0285962       2R   9864510-9875072      - | FBgn0285962     CG46338
##   FBgn0285963       3R 26915761-26945309      + | FBgn0285963     CG46339
##   FBgn0285970        X 21621541-21623750      - | FBgn0285970     CG32500
##   FBgn0285971       2L   8464488-8466694      + | FBgn0285971         prg
##                entrezid   gene_biotype seq_coord_system      symbol
##               <integer>    <character>        <integer> <character>
##   FBgn0000008      <NA> protein_coding             <NA>           a
##   FBgn0000014      <NA> protein_coding             <NA>       abd-A
##   FBgn0000015      <NA> protein_coding             <NA>       Abd-B
##   FBgn0000017      <NA> protein_coding             <NA>         Abl
##   FBgn0000018      <NA> protein_coding             <NA>         abo
##           ...       ...            ...              ...         ...
##   FBgn0285958      <NA> protein_coding             <NA>        Fuca
##   FBgn0285962      <NA> protein_coding             <NA>     CG46338
##   FBgn0285963      <NA> protein_coding             <NA>     CG46339
##   FBgn0285970      <NA> protein_coding             <NA>     CG32500
##   FBgn0285971      <NA> protein_coding             <NA>         prg
##   -------
##   seqinfo: 25 sequences from BDGP6 genome
```

# Add different identifiers

We would like to add support to easily map transcript or gene
identifiers from one annotation to another. This is just a prototype
function, but we show how we can easily add alternate IDs given that we
know the organism and the source of the transcriptome. (This function
currently only works for Gencode and Ensembl gene or transcript IDs
but could be extended to work for arbitrary sources.)


```r
library(org.Dm.eg.db)
gse <- addIds(gse, "REFSEQ", gene=TRUE)
```

```
## mapping to new IDs using 'org.Dm.eg.db' data package
```

```r
mcols(gse)
```

```
## DataFrame with 14026 rows and 7 columns
##           gene_id   gene_name  entrezid   gene_biotype seq_coord_system
##       <character> <character> <integer>    <character>        <integer>
## 1     FBgn0000008           a        NA protein_coding               NA
## 2     FBgn0000014       abd-A        NA protein_coding               NA
## 3     FBgn0000015       Abd-B        NA protein_coding               NA
## 4     FBgn0000017         Abl        NA protein_coding               NA
## 5     FBgn0000018         abo        NA protein_coding               NA
## ...           ...         ...       ...            ...              ...
## 14022 FBgn0285958        Fuca        NA protein_coding               NA
## 14023 FBgn0285962     CG46338        NA protein_coding               NA
## 14024 FBgn0285963     CG46339        NA protein_coding               NA
## 14025 FBgn0285970     CG32500        NA protein_coding               NA
## 14026 FBgn0285971         prg        NA protein_coding               NA
##            symbol       REFSEQ
##       <character>  <character>
## 1               a NM_001014543
## 2           abd-A NM_001170161
## 3           Abd-B NM_001275719
## 4             Abl NM_001104153
## 5             abo    NM_080045
## ...           ...          ...
## 14022        Fuca NM_001316434
## 14023     CG46338 NM_001273908
## 14024     CG46339 NM_001104469
## 14025     CG32500    NM_167766
## 14026         prg NM_001273324
```

# Run a differential expression analysis

The following code chunk demonstrates how to build a *DESeqDataSet*
and begin a differential expression analysis. Likely we would
simplify these steps with a convenience function either in *tximeta*
or in *DESeq2*.


```r
suppressPackageStartupMessages(library(DESeq2))
gse2 <- gse
assayNames(gse2)
```

```
## [1] "counts"    "abundance" "length"
```

```r
# rounding counts
assay(gse2) <- round(assay(gse2))
# rename the "length" assay to "avgTxLength"
# DESeq2 will then use this matrix for calculating bias offsets
assayNames(gse2)[3] <- "avgTxLength"
# here there is a single sample so we use ~1
dds <- DESeqDataSet(gse2, ~1)
```

```
## converting counts to integer mode
```

```
## Warning in DESeqDataSet(gse2, ~1): all genes have equal values for all
## samples. will not be able to perform differential analysis
```

```r
dds <- estimateSizeFactors(dds)
```

```
## using 'avgTxLength' from assays(dds), correcting for library size
```

```r
# ... and so on
```

# Metadata galore

The following information is attached to the *SummarizedExperiment* by
`tximeta`: 


```r
names(metadata(se))
```

```
## [1] "quantInfo"   "tximetaInfo" "txomeInfo"   "txdbInfo"
```

```r
str(metadata(se)$quantInfo)
```

```
## List of 24
##  $ salmon_version       : chr "0.10.2"
##  $ samp_type            : chr "none"
##  $ quant_errors         :List of 1
##   ..$ : list()
##  $ num_libraries        : int 1
##  $ library_types        : chr "ISR"
##  $ frag_dist_length     : int 1001
##  $ seq_bias_correct     : logi TRUE
##  $ gc_bias_correct      : logi TRUE
##  $ num_bias_bins        : int 4096
##  $ mapping_type         : chr "mapping"
##  $ num_targets          : int 30061
##  $ serialized_eq_classes: logi FALSE
##  $ length_classes       : int [1:5, 1] 1071 1736 2594 4068 71382
##  $ index_seq_hash       : chr "b41ea9ba9c81e2cad7cfa49e4bf9ee67dd297dc0b9ff40bdb1142699f00c8f7d"
##  $ index_name_hash      : chr "6aba201931d0fa4c6cebd3c1d7dd6350bf65cc1c968e88a308fe147f8a1c7083"
##  $ index_seq_hash512    : chr "365e58ceacde84989cb2bcc01e5b5c3320345ef23b23f1b49456f2ae429b5be5c418e233729d7a065e59645b5f26e25defbc1df40a601e6"| __truncated__
##  $ index_name_hash512   : chr "df81eababb8186637181132cf98221f9fa6fe77cb45ed771cc81bcfdda281cea5b89877ba233d81a06a1caa5aeedbe1b3161bd19c29bd6a"| __truncated__
##  $ num_bootstraps       : int 0
##  $ num_processed        : int 42422337
##  $ num_mapped           : int 29341160
##  $ percent_mapped       : num 69.2
##  $ call                 : chr "quant"
##  $ start_time           : chr "Fri Jul 13 08:45:38 2018"
##  $ end_time             : chr "Fri Jul 13 08:57:39 2018"
```

```r
str(metadata(se)$txomeInfo)
```

```
## List of 8
##  $ index         : chr "Drosophila_melanogaster.BDGP6.cdna.v92_salmon_0.10.2"
##  $ index_seq_hash: chr "b41ea9ba9c81e2cad7cfa49e4bf9ee67dd297dc0b9ff40bdb1142699f00c8f7d"
##  $ source        : chr "Ensembl"
##  $ organism      : chr "Drosophila melanogaster"
##  $ version       : chr "92"
##  $ genome        : chr "BDGP6"
##  $ fasta         :List of 1
##   ..$ : chr "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
##  $ gtf           : chr "/Users/love/Library/R/3.5/library/tximportData/extdata/salmon_dm/Drosophila_melanogaster.BDGP6.92.gtf.gz"
```

```r
str(metadata(se)$tximetaInfo)
```

```
## List of 2
##  $ version   :Classes 'package_version', 'numeric_version'  hidden list of 1
##   ..$ : int [1:3] 0 0 13
##  $ importTime: POSIXct[1:1], format: "2018-07-14 21:14:22"
```

```r
str(metadata(se)$txdbInfo)
```

```
##  Named chr [1:12] "EnsDb" "Ensembl Gene ID" "ensembldb" ...
##  - attr(*, "names")= chr [1:12] "Db type" "Type of Gene ID" "Supporting package" "Db created by" ...
```



# Next steps

### Basic functionality

* Switching `rowRanges` from transcript ranges to exons-by-transcript
  ranges list, or from gene ranges to exons-by-gene ranges list.
* As is already supported in `tximport`, also import inferential
  variance matrices (Gibbs samples or bootstrap samples)

### Facilitate plots and summaries
    
* Basic plots across samples: abundances, mapping rates, rich bias model parameters
* Time summaries: when quantified? when imported? I would love to
  know when the library was prepared and sequenced but this seems hopeless.

### Challenges

* Building out actual, sustainable plan for supporting as many
  organisms and sources as possible. We can define rules which
  determine where the FASTA and GTF files will be based on `source` and
  `version` (also here I ignored something like "type", e.g. CHR
  or ALL gene files from Gencode)
* Some support already for linked transcriptomes, see `linkedTxomes`
  vignette. Need to work more on combining multiple sources
  (potentially meta-transcriptomes from different organisms?), and
  also on how to approach de novo transcriptomes, and how to support
  reproducibility there.
* Facilitate functional annotation, either with vignettes/workflow or
  with additional functionality. E.g.: 
  housekeeping genes, arbitrary gene sets, genes expressed in GTEx tissues
* `liftOver` is clunky and doesn't integrate with
  *GenomeInfoDb*. It requires user input and there's a chance to
  mis-annotate. Ideally this should all be automated.

# Session info


```r
library(devtools)
session_info()
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.5.0 (2018-04-23)
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       Europe/Rome                 
##  date     2018-07-14
```

```
## Packages -----------------------------------------------------------------
```

```
##  package                           * version   date      
##  acepack                             1.4.1     2016-10-29
##  annotate                            1.58.0    2018-05-01
##  AnnotationDbi                     * 1.42.1    2018-05-08
##  AnnotationFilter                  * 1.4.0     2018-05-01
##  AnnotationHub                     * 2.12.0    2018-05-01
##  assertthat                          0.2.0     2017-04-11
##  backports                           1.1.2     2017-12-13
##  base                              * 3.5.0     2018-04-24
##  base64enc                           0.1-3     2015-07-28
##  bindr                               0.1.1     2018-03-13
##  bindrcpp                          * 0.2.2     2018-03-29
##  Biobase                           * 2.40.0    2018-05-01
##  BiocFileCache                     * 1.4.0     2018-05-01
##  BiocGenerics                      * 0.26.0    2018-05-01
##  BiocInstaller                     * 1.30.0    2018-05-04
##  BiocParallel                      * 1.14.1    2018-05-06
##  biomaRt                             2.36.1    2018-05-24
##  Biostrings                        * 2.48.0    2018-05-01
##  bit                                 1.1-14    2018-05-29
##  bit64                               0.9-7     2017-05-08
##  bitops                              1.0-6     2013-08-17
##  blob                                1.1.1     2018-03-25
##  checkmate                           1.8.5     2017-10-24
##  cluster                             2.0.7-1   2018-04-13
##  codetools                           0.2-15    2016-10-05
##  colorspace                          1.3-2     2016-12-14
##  commonmark                          1.5       2018-04-28
##  compiler                            3.5.0     2018-04-24
##  crayon                              1.3.4     2017-09-16
##  curl                                3.2       2018-03-28
##  data.table                          1.11.4    2018-05-27
##  datasets                          * 3.5.0     2018-04-24
##  DBI                                 1.0.0     2018-05-02
##  dbplyr                            * 1.2.1     2018-02-19
##  DelayedArray                      * 0.6.1     2018-06-15
##  DESeq2                            * 1.20.0    2018-05-01
##  devtools                          * 1.13.6    2018-06-27
##  digest                              0.6.15    2018-01-28
##  dplyr                               0.7.5     2018-05-19
##  ensembldb                         * 2.4.1     2018-05-07
##  evaluate                            0.10.1    2017-06-24
##  foreign                             0.8-70    2017-11-28
##  Formula                             1.2-3     2018-05-03
##  genefilter                          1.62.0    2018-05-01
##  geneplotter                         1.58.0    2018-05-01
##  GenomeInfoDb                      * 1.16.0    2018-05-01
##  GenomeInfoDbData                    1.1.0     2018-01-10
##  GenomicAlignments                   1.16.0    2018-05-01
##  GenomicFeatures                   * 1.32.0    2018-05-01
##  GenomicRanges                     * 1.32.3    2018-05-16
##  ggplot2                             3.0.0     2018-07-03
##  git2r                               0.21.0    2018-01-04
##  glue                                1.2.0     2017-10-29
##  GO.db                             * 3.5.0     2017-12-07
##  graph                               1.58.0    2018-05-01
##  graphics                          * 3.5.0     2018-04-24
##  grDevices                         * 3.5.0     2018-04-24
##  grid                                3.5.0     2018-04-24
##  gridExtra                           2.3       2017-09-09
##  gtable                              0.2.0     2016-02-26
##  Hmisc                               4.1-1     2018-01-03
##  hms                                 0.4.2     2018-03-10
##  Homo.sapiens                      * 1.3.1     2018-07-12
##  htmlTable                           1.12      2018-05-26
##  htmltools                           0.3.6     2017-04-28
##  htmlwidgets                         1.2       2018-04-19
##  httpuv                              1.4.3     2018-05-10
##  httr                                1.3.1     2017-08-20
##  interactiveDisplayBase              1.18.0    2018-05-01
##  IRanges                           * 2.14.10   2018-05-16
##  jsonlite                            1.5       2017-06-01
##  knitr                               1.20      2018-02-20
##  later                               0.7.2     2018-05-01
##  lattice                             0.20-35   2017-03-25
##  latticeExtra                        0.6-28    2016-02-09
##  lazyeval                            0.2.1     2017-10-29
##  locfit                              1.5-9.1   2013-04-20
##  magrittr                            1.5       2014-11-22
##  Matrix                              1.2-14    2018-04-13
##  matrixStats                       * 0.53.1    2018-02-11
##  memoise                             1.1.0     2017-04-21
##  methods                           * 3.5.0     2018-04-24
##  mime                                0.5       2016-07-07
##  munsell                             0.5.0     2018-06-12
##  nnet                                7.3-12    2016-02-02
##  org.Dm.eg.db                      * 3.6.0     2018-07-14
##  org.Hs.eg.db                      * 3.5.0     2017-12-07
##  OrganismDbi                       * 1.22.0    2018-05-01
##  parallel                          * 3.5.0     2018-04-24
##  pillar                              1.2.3     2018-05-25
##  pkgconfig                           2.0.1     2017-03-21
##  plyr                                1.8.4     2016-06-08
##  prettyunits                         1.0.2     2015-07-13
##  progress                            1.2.0     2018-06-14
##  promises                            1.0.1     2018-04-13
##  ProtGenerics                        1.12.0    2018-05-01
##  purrr                               0.2.5     2018-05-29
##  R6                                  2.2.2     2017-06-17
##  rappdirs                            0.3.1     2016-03-28
##  RBGL                                1.56.0    2018-05-01
##  RColorBrewer                        1.1-2     2014-12-07
##  Rcpp                                0.12.17   2018-05-18
##  RCurl                               1.95-4.10 2018-01-04
##  readr                               1.1.1     2017-05-16
##  rlang                               0.2.1     2018-05-30
##  rmarkdown                         * 1.9       2018-03-01
##  roxygen2                            6.0.1     2017-02-06
##  rpart                               4.1-13    2018-02-23
##  rprojroot                           1.3-2     2018-01-03
##  Rsamtools                           1.32.2    2018-07-03
##  RSQLite                             2.1.1     2018-05-06
##  rstudioapi                          0.7       2017-09-07
##  rtracklayer                       * 1.40.3    2018-06-02
##  S4Vectors                         * 0.18.3    2018-06-08
##  scales                              0.5.0     2017-08-24
##  shiny                               1.0.5     2017-08-23
##  splines                             3.5.0     2018-04-24
##  stats                             * 3.5.0     2018-04-24
##  stats4                            * 3.5.0     2018-04-24
##  stringi                             1.2.3     2018-06-12
##  stringr                             1.3.1     2018-05-10
##  SummarizedExperiment              * 1.10.1    2018-05-11
##  survival                            2.42-3    2018-04-16
##  testthat                          * 2.0.0     2017-12-13
##  tibble                              1.4.2     2018-01-22
##  tidyselect                          0.2.4     2018-02-26
##  tools                               3.5.0     2018-04-24
##  TxDb.Hsapiens.UCSC.hg19.knownGene * 3.2.2     2018-05-04
##  tximeta                           * 0.0.13    2018-07-14
##  tximport                          * 1.8.0     2018-05-01
##  utils                             * 3.5.0     2018-04-24
##  withr                               2.1.2     2018-03-15
##  XML                                 3.98-1.11 2018-04-16
##  xml2                                1.2.0     2018-01-24
##  xtable                              1.8-2     2016-02-05
##  XVector                           * 0.20.0    2018-05-01
##  yaml                                2.1.19    2018-05-01
##  zlibbioc                            1.26.0    2018-05-01
##  source                     
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  cran (@1.1.2)              
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  cran (@0.6.15)             
##  cran (@0.7.5)              
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  Bioconductor               
##  local                      
##  local                      
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  Bioconductor               
##  Bioconductor               
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  cran (@0.2.5)              
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  cran (@0.12.17)            
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  cran (@0.2.1)              
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  cran (@1.3-2)              
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  cran (@1.40.3)             
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  local                      
##  local                      
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  local                      
##  Bioconductor               
##  local (mikelove/tximeta@NA)
##  Bioconductor               
##  local                      
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  CRAN (R 3.5.0)             
##  Bioconductor               
##  CRAN (R 3.5.0)             
##  Bioconductor
```
