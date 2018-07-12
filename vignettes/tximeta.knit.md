---
title: "tximeta: Import transcript quantification with automagic generation of metadata"
author: "Michael Love, Rob Patro"
date: "07/12/2018"
output: 
  html_document:
    highlight: tango
abstract: >
  `tximeta` performs numerous annotation and metadata gathering tasks on
  behalf of users during the import of transcript quantifications from
  Salmon or Sailfish into R/Bioconductor. The goal is to provide
  something similar to the experience of `GEOquery`, which downloaded
  microarray expression data from NCBI GEO and simultaneously brought
  along associated pieces of metadata. Doing this automatically helps to
  prevent costly bioinformatic errors. 
---

# This package is in beta 

[See README](https://github.com/mikelove/tximeta/blob/master/README.md)



# Analysis starts with sample table

The first step using `tximeta` is to read in the sample table, which
will become the *column data*, `colData`, of the
*SummarizedExperiment*. This table should contain all the information
we need to identify the quant directories. Here we will use the
*Salmon* quantification files in the *tximportData* package to
demonstrate the usage of `tximeta`. We do not have a sample table, so
we construct one in R, but it is recommended to keep a sample table as
a CSV or TSV file while working on a project.


```r
dir <- system.file("extdata/salmon", package="tximportData")
ids <- list.files(dir)
# here gzipped, normally these are not
files <- file.path(dir, ids, "quant.sf.gz") 
file.exists(files)
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE
```

```r
coldata <- data.frame(files, names=ids, population="TSI", stringsAsFactors=FALSE)
coldata
```

```
##                                                                                 files
## 1 /Users/love/Library/R/3.5/library/tximportData/extdata/salmon/ERR188021/quant.sf.gz
## 2 /Users/love/Library/R/3.5/library/tximportData/extdata/salmon/ERR188088/quant.sf.gz
## 3 /Users/love/Library/R/3.5/library/tximportData/extdata/salmon/ERR188288/quant.sf.gz
## 4 /Users/love/Library/R/3.5/library/tximportData/extdata/salmon/ERR188297/quant.sf.gz
## 5 /Users/love/Library/R/3.5/library/tximportData/extdata/salmon/ERR188329/quant.sf.gz
## 6 /Users/love/Library/R/3.5/library/tximportData/extdata/salmon/ERR188356/quant.sf.gz
##       names population
## 1 ERR188021        TSI
## 2 ERR188088        TSI
## 3 ERR188288        TSI
## 4 ERR188297        TSI
## 5 ERR188329        TSI
## 6 ERR188356        TSI
```

`tximeta` expects at least two columns in `coldata`: 

1. `files` - a pointer to the `quant.sf` files
2. `names` - the unique names that should be used to identify samples

# Running tximeta from a sample table

(Note: first do a `devtools::load_all()` then the following should work.)


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
## 1 2 3 4 5 6 
## found matching transcriptome:
## [ Gencode - Homo sapiens - version 27 ]
## loading existing TxDb created: 2018-07-12 13:44:25
## generating transcript ranges
## fetching genome info
```

# What happened? 

`tximeta` recognized the signature of the transcriptome
that the files were quantified against, it accessed the remote GTF
file of the transcriptome source, found and attached the transcript
genomic ranges, and added the appropriate transcriptome and genome metadata.
The remote GTF is only accessed once. If `tximeta` recognizes that it
has seen this index before, it will simply use a cached version of the
transcript metadata (this uses *BiocFileCache* and the specifics of
the cache location may change as `tximeta` develops).

We plan to create and maintain a large table of signatures for as many
sources, organisms, versions of transcriptomes as possible. We are
also developing support for "derived transcriptomes", where one or
more sources for transcript sequences have been merged or
filtered. See the `derivedTxome` vignette in this package for a
demonstration. 

# Examining SummarizedExperiment output

We, of course, have our coldata from before. Note that we've removed `files`.


```r
colData(se)
```

```
## DataFrame with 6 rows and 2 columns
##                 names  population
##           <character> <character>
## ERR188021   ERR188021         TSI
## ERR188088   ERR188088         TSI
## ERR188288   ERR188288         TSI
## ERR188297   ERR188297         TSI
## ERR188329   ERR188329         TSI
## ERR188356   ERR188356         TSI
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

Thanks to `tximeta` we have automagically imported the correct ranges
for the transcripts. 


```r
rowRanges(se)
```

```
## GRanges object with 200401 ranges and 2 metadata columns:
##                     seqnames      ranges strand |     tx_id
##                        <Rle>   <IRanges>  <Rle> | <integer>
##   ENST00000456328.2     chr1 11869-14409      + |         1
##   ENST00000450305.2     chr1 12010-13670      + |         2
##   ENST00000488147.1     chr1 14404-29570      - |      9081
##   ENST00000619216.1     chr1 17369-17436      - |      9082
##   ENST00000473358.1     chr1 29554-31097      + |         3
##                 ...      ...         ...    ... .       ...
##   ENST00000361681.2     chrM 14149-14673      - |    200399
##   ENST00000387459.1     chrM 14674-14742      - |    200400
##   ENST00000361789.2     chrM 14747-15887      + |    200391
##   ENST00000387460.2     chrM 15888-15953      + |    200392
##   ENST00000387461.2     chrM 15956-16023      - |    200401
##                               tx_name
##                           <character>
##   ENST00000456328.2 ENST00000456328.2
##   ENST00000450305.2 ENST00000450305.2
##   ENST00000488147.1 ENST00000488147.1
##   ENST00000619216.1 ENST00000619216.1
##   ENST00000473358.1 ENST00000473358.1
##                 ...               ...
##   ENST00000361681.2 ENST00000361681.2
##   ENST00000387459.1 ENST00000387459.1
##   ENST00000361789.2 ENST00000361789.2
##   ENST00000387460.2 ENST00000387460.2
##   ENST00000387461.2 ENST00000387461.2
##   -------
##   seqinfo: 25 sequences (1 circular) from hg38 genome
```

We have appropriate genome information, which prevents us from making 
bioinformatic mistakes:


```r
seqinfo(se)
```

```
## Seqinfo object with 25 sequences (1 circular) from hg38 genome:
##   seqnames seqlengths isCircular genome
##   chr1      248956422      FALSE   hg38
##   chr2      242193529      FALSE   hg38
##   chr3      198295559      FALSE   hg38
##   chr4      190214555      FALSE   hg38
##   chr5      181538259      FALSE   hg38
##   ...             ...        ...    ...
##   chr21      46709983      FALSE   hg38
##   chr22      50818468      FALSE   hg38
##   chrX      156040895      FALSE   hg38
##   chrY       57227415      FALSE   hg38
##   chrM          16569       TRUE   hg38
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
## loading existing TxDb created: 2018-07-12 13:44:25
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
## GRanges object with 58288 ranges and 1 metadata column:
##                      seqnames              ranges strand |
##                         <Rle>           <IRanges>  <Rle> |
##   ENSG00000000003.14     chrX 100627109-100639991      - |
##    ENSG00000000005.5     chrX 100584802-100599885      + |
##   ENSG00000000419.12    chr20   50934867-50958555      - |
##   ENSG00000000457.13     chr1 169849631-169894267      - |
##   ENSG00000000460.16     chr1 169662007-169854080      + |
##                  ...      ...                 ...    ... .
##    ENSG00000284744.1     chr1     6767954-6770038      + |
##    ENSG00000284745.1     chr1     2960658-2968707      - |
##    ENSG00000284746.1     chr8   12601158-12601376      - |
##    ENSG00000284747.1     chr1     7991134-8005312      + |
##    ENSG00000284748.1     chr1   37596126-37607336      + |
##                                 gene_id
##                             <character>
##   ENSG00000000003.14 ENSG00000000003.14
##    ENSG00000000005.5  ENSG00000000005.5
##   ENSG00000000419.12 ENSG00000000419.12
##   ENSG00000000457.13 ENSG00000000457.13
##   ENSG00000000460.16 ENSG00000000460.16
##                  ...                ...
##    ENSG00000284744.1  ENSG00000284744.1
##    ENSG00000284745.1  ENSG00000284745.1
##    ENSG00000284746.1  ENSG00000284746.1
##    ENSG00000284747.1  ENSG00000284747.1
##    ENSG00000284748.1  ENSG00000284748.1
##   -------
##   seqinfo: 25 sequences (1 circular) from an unspecified genome; no seqlengths
```

# Add different identifiers

We would like to add support to easily map transcript or gene
identifiers from one annotation to another. This is just a prototype
function, but we show how we can easily add alternate IDs given that we
know the organism and the source of the transcriptome. (This function
currently only works for Gencode and Ensembl gene or transcript IDs
for human, but could be extended to work for arbitrary sources.)


```r
library(Homo.sapiens)
```

```
## Loading required package: OrganismDbi
```

```
## Loading required package: GO.db
```

```
## 
```

```
## Loading required package: org.Hs.eg.db
```

```
## 
```

```
## Loading required package: TxDb.Hsapiens.UCSC.hg19.knownGene
```

```r
gse <- addIds(gse, "SYMBOL")
```

```
## mapping to new IDs using 'Homo.sapiens' data package
```

```r
mcols(gse)[1:10,]
```

```
## DataFrame with 10 rows and 2 columns
##               gene_id      SYMBOL
##           <character> <character>
## 1  ENSG00000000003.14      TSPAN6
## 2   ENSG00000000005.5        TNMD
## 3  ENSG00000000419.12        DPM1
## 4  ENSG00000000457.13       SCYL3
## 5  ENSG00000000460.16    C1orf112
## 6  ENSG00000000938.12         FGR
## 7  ENSG00000000971.15         CFH
## 8  ENSG00000001036.13       FUCA2
## 9  ENSG00000001084.10        GCLC
## 10 ENSG00000001167.14        NFYA
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
# make a DESeqDataSet, here there is no factor
# to divide the samples so we use ~1
dds <- DESeqDataSet(gse2, ~1)
```

```
## converting counts to integer mode
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

# Find nearest transcripts to a ChIP-seq peak

Suppose we want to find overlap of the expression with binding sites
of a transcription factor:


```r
library(AnnotationHub)
```

```
## 
## Attaching package: 'AnnotationHub'
```

```
## The following object is masked from 'package:Biobase':
## 
##     cache
```

```r
ah <- AnnotationHub()
```

```
## snapshotDate(): 2018-04-30
```

```r
chip <- query(ah, c("GM12878", "MEF2A", "narrowPeak"))[[1]]
```

```
## require("rtracklayer")
```

```
## downloading 0 resources
```

```
## loading from cache 
##     '/Users/love//.AnnotationHub/28040'
```

First try, let's find the nearest transcript to a given ChIP-seq peak:


```r
nearest(chip[1], se)
```

```
## Error in mergeNamedAtomicVectors(genome(x), genome(y), what = c("sequence", : sequences chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY, chrM have incompatible genomes:
##   - in 'x': hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19, hg19
##   - in 'y': hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38, hg38
```

We get an <font color="red"><b>ERROR</b></font>: all chromosomes have
incompatibile genomes! Good! 
The point of `tximeta` is to reduce these kind of simple bioinformatic
mistakes that can add weeks or months of dead-end results to large
genomics projects. 
We can use liftover chains to get us going in the right direction,
comparing hg38 to hg38 (this code chunk un-evaluated).


```r
url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
file <- "hg19ToHg38.over.chain.gz"
if (!file.exists(file)) download.file(url, file, method="wget")
system(paste("gunzip ",file))
```


```r
chainfile <- system.file("extdata/hg19ToHg38.over.chain", package="tximeta")
```

We move our ChIP-seq data to hg38:


```r
chip.lift <- liftOverHelper(chip, chainfile=chainfile, to="hg38")
```

Now we can find the nearest transcript to a given ChIP-seq peak:


```r
nearest(chip.lift[1], se)
```

```
## [1] 62266
```

```r
assay(se)[nearest(chip.lift[1], se),,drop=FALSE]
```

```
##                   ERR188021 ERR188088 ERR188288 ERR188297 ERR188329
## ENST00000379350.5        41   15.0285   21.8972   21.2163   4.02264
##                   ERR188356
## ENST00000379350.5   20.0937
```

Or we can take a slice of the transcriptome data that is within 10 kb
of a given ChIP-seq peak:


```r
# which rows of SE in this window?
which(overlapsAny(se, chip.lift[1] + 1e4))
```

```
## [1] 62265 62266 62276 62280
```

Perhaps even more exciting, we can now automate functional annotation of
transcriptome data using Bioconductor's annotation suite.

# Metadata galore


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
## List of 21
##  $ salmon_version       : chr [1:6] "0.8.2" "0.8.2" "0.8.2" "0.8.2" ...
##  $ samp_type            : chr [1:6] "none" "none" "none" "none" ...
##  $ num_libraries        : int [1:6] 1 1 1 1 1 1
##  $ library_types        : chr [1:6] "IU" "IU" "IU" "IU" ...
##  $ frag_dist_length     : int [1:6] 1001 1001 1001 1001 1001 1001
##  $ seq_bias_correct     : logi [1:6] FALSE FALSE FALSE FALSE FALSE FALSE
##  $ gc_bias_correct      : logi [1:6] TRUE TRUE TRUE TRUE TRUE TRUE
##  $ num_bias_bins        : int [1:6] 4096 4096 4096 4096 4096 4096
##  $ mapping_type         : chr [1:6] "mapping" "mapping" "mapping" "mapping" ...
##  $ num_targets          : int [1:6] 200401 200401 200401 200401 200401 200401
##  $ serialized_eq_classes: logi [1:6] FALSE FALSE FALSE FALSE FALSE FALSE
##  $ length_classes       : int [1:5, 1:6] 513 655 1015 2261 205012 513 655 1015 2261 205012 ...
##  $ index_seq_hash       : chr "0acbb6919651beb0182512f492e9936c43dd7c0d189df41c402fbd4fc556f460"
##  $ index_name_hash      : chr "92aad1c26de97ac5d9de6cbbc66424ac8bfc1b1ade3be841d586e57076ac5b58"
##  $ num_bootstraps       : int [1:6] 0 0 0 0 0 0
##  $ num_processed        : int [1:6] 32507828 26194108 27045305 25288232 37938491 26304702
##  $ num_mapped           : int [1:6] 28694462 23721986 23854225 22830167 34481916 22602902
##  $ percent_mapped       : num [1:6] 88.3 90.6 88.2 90.3 90.9 ...
##  $ call                 : chr [1:6] "quant" "quant" "quant" "quant" ...
##  $ start_time           : chr [1:6] "Wed Dec  6 16:56:54 2017" "Wed Dec  6 16:16:00 2017" "Wed Dec  6 16:41:41 2017" "Wed Dec  6 16:03:00 2017" ...
##  $ end_time             : chr [1:6] "Wed Dec  6 17:10:44 2017" "Wed Dec  6 16:27:50 2017" "Wed Dec  6 16:56:53 2017" "Wed Dec  6 16:15:59 2017" ...
```

```r
str(metadata(se)$txomeInfo)
```

```
## List of 7
##  $ index_seq_hash: chr "0acbb6919651beb0182512f492e9936c43dd7c0d189df41c402fbd4fc556f460"
##  $ source        : chr "Gencode"
##  $ organism      : chr "Homo sapiens"
##  $ version       : chr "27"
##  $ genome        : chr "GRCh38"
##  $ fasta         : chr "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz"
##  $ gtf           : chr "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz"
```

```r
str(metadata(se)$tximetaInfo)
```

```
## List of 2
##  $ version   :Classes 'package_version', 'numeric_version'  hidden list of 1
##   ..$ : int [1:3] 0 0 10
##  $ importTime: POSIXct[1:1], format: "2018-07-12 16:30:57"
```

```r
str(metadata(se)$txdbInfo)
```

```
##  Named chr [1:15] "TxDb" "GenomicFeatures" ...
##  - attr(*, "names")= chr [1:15] "Db type" "Supporting package" "Data source" "Organism" ...
```

# Next steps

### Basic functionality

* Switching `rowRanges` from transcript ranges to exons-by-transcript
  ranges list, or from gene ranges to exons-by-gene ranges list.
* As is already supported in the `tximport` release, also import inferential
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
* Some support already for derived transcriptomes, see `derivedTxomes`
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
##  date     2018-07-12
```

```
## Packages -----------------------------------------------------------------
```

```
##  package                           * version   date       source         
##  acepack                             1.4.1     2016-10-29 CRAN (R 3.5.0) 
##  annotate                            1.58.0    2018-05-01 Bioconductor   
##  AnnotationDbi                     * 1.42.1    2018-05-08 Bioconductor   
##  AnnotationFilter                  * 1.4.0     2018-05-01 Bioconductor   
##  AnnotationHub                     * 2.12.0    2018-05-01 Bioconductor   
##  assertthat                          0.2.0     2017-04-11 CRAN (R 3.5.0) 
##  backports                           1.1.2     2017-12-13 cran (@1.1.2)  
##  base                              * 3.5.0     2018-04-24 local          
##  base64enc                           0.1-3     2015-07-28 CRAN (R 3.5.0) 
##  bindr                               0.1.1     2018-03-13 CRAN (R 3.5.0) 
##  bindrcpp                          * 0.2.2     2018-03-29 CRAN (R 3.5.0) 
##  Biobase                           * 2.40.0    2018-05-01 Bioconductor   
##  BiocFileCache                     * 1.4.0     2018-05-01 Bioconductor   
##  BiocGenerics                      * 0.26.0    2018-05-01 Bioconductor   
##  BiocInstaller                     * 1.30.0    2018-05-04 Bioconductor   
##  BiocParallel                      * 1.14.1    2018-05-06 Bioconductor   
##  biomaRt                             2.36.1    2018-05-24 Bioconductor   
##  Biostrings                        * 2.48.0    2018-05-01 Bioconductor   
##  bit                                 1.1-14    2018-05-29 CRAN (R 3.5.0) 
##  bit64                               0.9-7     2017-05-08 CRAN (R 3.5.0) 
##  bitops                              1.0-6     2013-08-17 CRAN (R 3.5.0) 
##  blob                                1.1.1     2018-03-25 CRAN (R 3.5.0) 
##  BSgenome                            1.48.0    2018-05-01 Bioconductor   
##  checkmate                           1.8.5     2017-10-24 CRAN (R 3.5.0) 
##  cli                                 1.0.0     2017-11-05 CRAN (R 3.5.0) 
##  cluster                             2.0.7-1   2018-04-13 CRAN (R 3.5.0) 
##  codetools                           0.2-15    2016-10-05 CRAN (R 3.5.0) 
##  colorspace                          1.3-2     2016-12-14 CRAN (R 3.5.0) 
##  commonmark                          1.5       2018-04-28 CRAN (R 3.5.0) 
##  compiler                            3.5.0     2018-04-24 local          
##  crayon                              1.3.4     2017-09-16 CRAN (R 3.5.0) 
##  curl                                3.2       2018-03-28 CRAN (R 3.5.0) 
##  data.table                          1.11.4    2018-05-27 CRAN (R 3.5.0) 
##  datasets                          * 3.5.0     2018-04-24 local          
##  DBI                                 1.0.0     2018-05-02 CRAN (R 3.5.0) 
##  dbplyr                            * 1.2.1     2018-02-19 CRAN (R 3.5.0) 
##  DelayedArray                      * 0.6.1     2018-06-15 Bioconductor   
##  desc                                1.2.0     2018-05-01 CRAN (R 3.5.0) 
##  DESeq2                            * 1.20.0    2018-05-01 Bioconductor   
##  devtools                          * 1.13.6    2018-06-27 CRAN (R 3.5.0) 
##  digest                              0.6.15    2018-01-28 cran (@0.6.15) 
##  dplyr                               0.7.5     2018-05-19 cran (@0.7.5)  
##  ensembldb                         * 2.4.1     2018-05-07 Bioconductor   
##  evaluate                            0.10.1    2017-06-24 CRAN (R 3.5.0) 
##  foreign                             0.8-70    2017-11-28 CRAN (R 3.5.0) 
##  Formula                             1.2-3     2018-05-03 CRAN (R 3.5.0) 
##  genefilter                          1.62.0    2018-05-01 Bioconductor   
##  geneplotter                         1.58.0    2018-05-01 Bioconductor   
##  GenomeInfoDb                      * 1.16.0    2018-05-01 Bioconductor   
##  GenomeInfoDbData                    1.1.0     2018-01-10 Bioconductor   
##  GenomicAlignments                   1.16.0    2018-05-01 Bioconductor   
##  GenomicFeatures                   * 1.32.0    2018-05-01 Bioconductor   
##  GenomicRanges                     * 1.32.3    2018-05-16 Bioconductor   
##  ggplot2                             3.0.0     2018-07-03 CRAN (R 3.5.0) 
##  glue                                1.2.0     2017-10-29 CRAN (R 3.5.0) 
##  GO.db                             * 3.5.0     2017-12-07 Bioconductor   
##  graph                               1.58.0    2018-05-01 Bioconductor   
##  graphics                          * 3.5.0     2018-04-24 local          
##  grDevices                         * 3.5.0     2018-04-24 local          
##  grid                                3.5.0     2018-04-24 local          
##  gridExtra                           2.3       2017-09-09 CRAN (R 3.5.0) 
##  gtable                              0.2.0     2016-02-26 CRAN (R 3.5.0) 
##  Hmisc                               4.1-1     2018-01-03 CRAN (R 3.5.0) 
##  hms                                 0.4.2     2018-03-10 CRAN (R 3.5.0) 
##  Homo.sapiens                      * 1.3.1     2018-07-12 Bioconductor   
##  htmlTable                           1.12      2018-05-26 CRAN (R 3.5.0) 
##  htmltools                           0.3.6     2017-04-28 CRAN (R 3.5.0) 
##  htmlwidgets                         1.2       2018-04-19 CRAN (R 3.5.0) 
##  httpuv                              1.4.3     2018-05-10 CRAN (R 3.5.0) 
##  httr                                1.3.1     2017-08-20 CRAN (R 3.5.0) 
##  interactiveDisplayBase              1.18.0    2018-05-01 Bioconductor   
##  IRanges                           * 2.14.10   2018-05-16 Bioconductor   
##  jsonlite                            1.5       2017-06-01 CRAN (R 3.5.0) 
##  knitr                               1.20      2018-02-20 CRAN (R 3.5.0) 
##  later                               0.7.2     2018-05-01 CRAN (R 3.5.0) 
##  lattice                             0.20-35   2017-03-25 CRAN (R 3.5.0) 
##  latticeExtra                        0.6-28    2016-02-09 CRAN (R 3.5.0) 
##  lazyeval                            0.2.1     2017-10-29 CRAN (R 3.5.0) 
##  locfit                              1.5-9.1   2013-04-20 CRAN (R 3.5.0) 
##  magrittr                            1.5       2014-11-22 CRAN (R 3.5.0) 
##  Matrix                              1.2-14    2018-04-13 CRAN (R 3.5.0) 
##  matrixStats                       * 0.53.1    2018-02-11 CRAN (R 3.5.0) 
##  memoise                             1.1.0     2017-04-21 CRAN (R 3.5.0) 
##  methods                           * 3.5.0     2018-04-24 local          
##  mime                                0.5       2016-07-07 CRAN (R 3.5.0) 
##  munsell                             0.5.0     2018-06-12 CRAN (R 3.5.0) 
##  nnet                                7.3-12    2016-02-02 CRAN (R 3.5.0) 
##  org.Hs.eg.db                      * 3.5.0     2017-12-07 Bioconductor   
##  OrganismDbi                       * 1.22.0    2018-05-01 Bioconductor   
##  parallel                          * 3.5.0     2018-04-24 local          
##  pillar                              1.2.3     2018-05-25 CRAN (R 3.5.0) 
##  pkgconfig                           2.0.1     2017-03-21 CRAN (R 3.5.0) 
##  plyr                                1.8.4     2016-06-08 CRAN (R 3.5.0) 
##  prettyunits                         1.0.2     2015-07-13 CRAN (R 3.5.0) 
##  progress                            1.2.0     2018-06-14 CRAN (R 3.5.0) 
##  promises                            1.0.1     2018-04-13 CRAN (R 3.5.0) 
##  ProtGenerics                        1.12.0    2018-05-01 Bioconductor   
##  purrr                               0.2.5     2018-05-29 cran (@0.2.5)  
##  R6                                  2.2.2     2017-06-17 CRAN (R 3.5.0) 
##  rappdirs                            0.3.1     2016-03-28 CRAN (R 3.5.0) 
##  RBGL                                1.56.0    2018-05-01 Bioconductor   
##  RColorBrewer                        1.1-2     2014-12-07 CRAN (R 3.5.0) 
##  Rcpp                                0.12.17   2018-05-18 cran (@0.12.17)
##  RCurl                               1.95-4.10 2018-01-04 CRAN (R 3.5.0) 
##  readr                               1.1.1     2017-05-16 CRAN (R 3.5.0) 
##  rlang                               0.2.1     2018-05-30 cran (@0.2.1)  
##  rmarkdown                         * 1.9       2018-03-01 CRAN (R 3.5.0) 
##  roxygen2                            6.0.1     2017-02-06 CRAN (R 3.5.0) 
##  rpart                               4.1-13    2018-02-23 CRAN (R 3.5.0) 
##  rprojroot                           1.3-2     2018-01-03 cran (@1.3-2)  
##  Rsamtools                           1.32.2    2018-07-03 Bioconductor   
##  RSQLite                             2.1.1     2018-05-06 CRAN (R 3.5.0) 
##  rstudioapi                          0.7       2017-09-07 CRAN (R 3.5.0) 
##  rtracklayer                       * 1.40.2    2018-05-08 Bioconductor   
##  S4Vectors                         * 0.18.3    2018-06-08 Bioconductor   
##  scales                              0.5.0     2017-08-24 CRAN (R 3.5.0) 
##  shiny                               1.0.5     2017-08-23 CRAN (R 3.5.0) 
##  splines                             3.5.0     2018-04-24 local          
##  stats                             * 3.5.0     2018-04-24 local          
##  stats4                            * 3.5.0     2018-04-24 local          
##  stringi                             1.2.3     2018-06-12 CRAN (R 3.5.0) 
##  stringr                             1.3.1     2018-05-10 CRAN (R 3.5.0) 
##  SummarizedExperiment              * 1.10.1    2018-05-11 Bioconductor   
##  survival                            2.42-3    2018-04-16 CRAN (R 3.5.0) 
##  testthat                          * 2.0.0     2017-12-13 CRAN (R 3.5.0) 
##  tibble                              1.4.2     2018-01-22 CRAN (R 3.5.0) 
##  tidyselect                          0.2.4     2018-02-26 CRAN (R 3.5.0) 
##  tools                               3.5.0     2018-04-24 local          
##  TxDb.Hsapiens.UCSC.hg19.knownGene * 3.2.2     2018-05-04 Bioconductor   
##  tximeta                           * 0.0.10    <NA>       Bioconductor   
##  tximport                          * 1.8.0     2018-05-01 Bioconductor   
##  utf8                                1.1.4     2018-05-24 CRAN (R 3.5.0) 
##  utils                             * 3.5.0     2018-04-24 local          
##  withr                               2.1.2     2018-03-15 CRAN (R 3.5.0) 
##  XML                                 3.98-1.11 2018-04-16 CRAN (R 3.5.0) 
##  xml2                                1.2.0     2018-01-24 CRAN (R 3.5.0) 
##  xtable                              1.8-2     2016-02-05 CRAN (R 3.5.0) 
##  XVector                           * 0.20.0    2018-05-01 Bioconductor   
##  yaml                                2.1.19    2018-05-01 CRAN (R 3.5.0) 
##  zlibbioc                            1.26.0    2018-05-01 Bioconductor
```
