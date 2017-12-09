---
title: "tximeta: Import transcript quantification with automagic generation of metadata"
author: "Michael Love, Rob Patro"
date: "12/06/2017"
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

# Setup

First, to try out `tximeta` you'll need the example data, which is
contained in this GitHub repo. Here we download the whole repo as a
ZIP file.




```r
library(here)
dest <- here("extdata","asthma.zip")
download.file("https://github.com/mikelove/asthma/archive/master.zip", dest, method="wget")
unzip(dest, exdir=here("extdata"))
```

# Analysis starts with sample table

The first step using `tximeta` is to read in the sample table, which
will become the *column data*, `colData`, of the
*SummarizedExperiment*. This table should contain all the information
we need to identify the quant directories. In this case, we will use
the `run` ID from SRA.


```r
library(here)
library(readr)
coldata <- read_tsv(here("extdata","coldata.tsv"))
```

```
## Parsed with column specification:
## cols(
##   run = col_character(),
##   id = col_character(),
##   condition = col_character(),
##   treatment = col_character()
## )
```

```r
coldata
```

```
## # A tibble: 24 x 4
##           run    id condition treatment
##         <chr> <chr>     <chr>     <chr>
##  1 SRR1565944   s21       non   Vehicle
##  2 SRR1565945   s23       non   Vehicle
##  3 SRR1565946   s20       non   Vehicle
##  4 SRR1565947   s22       non   Vehicle
##  5 SRR1565948   s11       non   Vehicle
##  6 SRR1565949   s30       non   Vehicle
##  7 SRR1565943   s21       non     HRV16
##  8 SRR1565939   s23       non     HRV16
##  9 SRR1565940   s20       non     HRV16
## 10 SRR1565942   s22       non     HRV16
## # ... with 14 more rows
```

`tximeta` expects at least two columns in `coldata`: 

1. `files` - a pointer to the `quant.sf` files
2. `names` - the unique names that should be used to identify samples

We use `coldata$run` to build these two columns:


```r
coldata$files <- here("extdata","asthma-master","data","quants",
                      coldata$run,"quant.sf.gz")
coldata$names <- coldata$run
all(file.exists(coldata$files))
```

```
## [1] TRUE
```

```r
head(coldata$names)
```

```
## [1] "SRR1565944" "SRR1565945" "SRR1565946" "SRR1565947" "SRR1565948"
## [6] "SRR1565949"
```

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
## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 
## found matching transcriptome:
## [ Gencode - Homo sapiens - version 26 ]
## loading existing TxDb created: 2017-11-28 21:36:55
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
## DataFrame with 24 rows and 5 columns
##                    run          id   condition   treatment       names
##            <character> <character> <character> <character> <character>
## SRR1565944  SRR1565944         s21         non     Vehicle  SRR1565944
## SRR1565945  SRR1565945         s23         non     Vehicle  SRR1565945
## SRR1565946  SRR1565946         s20         non     Vehicle  SRR1565946
## SRR1565947  SRR1565947         s22         non     Vehicle  SRR1565947
## SRR1565948  SRR1565948         s11         non     Vehicle  SRR1565948
## ...                ...         ...         ...         ...         ...
## SRR1565926  SRR1565926         s12        asth       HRV16  SRR1565926
## SRR1565931  SRR1565931         s26        asth       HRV16  SRR1565931
## SRR1565930  SRR1565930         s15        asth       HRV16  SRR1565930
## SRR1565929  SRR1565929         s18        asth       HRV16  SRR1565929
## SRR1565927  SRR1565927         s14        asth       HRV16  SRR1565927
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
## GRanges object with 199324 ranges and 2 metadata columns:
##                     seqnames         ranges strand |     tx_id
##                        <Rle>      <IRanges>  <Rle> | <integer>
##   ENST00000456328.2     chr1 [11869, 14409]      + |         1
##   ENST00000450305.2     chr1 [12010, 13670]      + |         2
##   ENST00000488147.1     chr1 [14404, 29570]      - |      8972
##   ENST00000619216.1     chr1 [17369, 17436]      - |      8973
##   ENST00000473358.1     chr1 [29554, 31097]      + |         3
##                 ...      ...            ...    ... .       ...
##   ENST00000361681.2     chrM [14149, 14673]      - |    199322
##   ENST00000387459.1     chrM [14674, 14742]      - |    199323
##   ENST00000361789.2     chrM [14747, 15887]      + |    199314
##   ENST00000387460.2     chrM [15888, 15953]      + |    199315
##   ENST00000387461.2     chrM [15956, 16023]      - |    199324
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
## loading existing TxDb created: 2017-11-28 21:36:55
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
## GRanges object with 58219 ranges and 1 metadata column:
##                      seqnames                 ranges strand |
##                         <Rle>              <IRanges>  <Rle> |
##   ENSG00000000003.14     chrX [100627109, 100639991]      - |
##    ENSG00000000005.5     chrX [100584802, 100599885]      + |
##   ENSG00000000419.12    chr20 [ 50934867,  50958555]      - |
##   ENSG00000000457.13     chr1 [169849631, 169894267]      - |
##   ENSG00000000460.16     chr1 [169662007, 169854080]      + |
##                  ...      ...                    ...    ... .
##    ENSG00000284592.1     chr1 [157203604, 157205062]      - |
##    ENSG00000284594.1    chr11 [  1880045,   1880147]      + |
##    ENSG00000284595.1    chr17 [ 75498548,  75498628]      + |
##    ENSG00000284596.1     chr7 [102471469, 102471531]      + |
##    ENSG00000284600.1     chr2 [  1795525,   1811526]      + |
##                                 gene_id
##                             <character>
##   ENSG00000000003.14 ENSG00000000003.14
##    ENSG00000000005.5  ENSG00000000005.5
##   ENSG00000000419.12 ENSG00000000419.12
##   ENSG00000000457.13 ENSG00000000457.13
##   ENSG00000000460.16 ENSG00000000460.16
##                  ...                ...
##    ENSG00000284592.1  ENSG00000284592.1
##    ENSG00000284594.1  ENSG00000284594.1
##    ENSG00000284595.1  ENSG00000284595.1
##    ENSG00000284596.1  ENSG00000284596.1
##    ENSG00000284600.1  ENSG00000284600.1
##   -------
##   seqinfo: 25 sequences (1 circular) from an unspecified genome; no seqlengths
```

# Add different identifiers

We would like to add support to easily map transcript or gene
identifiers from one annotation to another. This is just a prototype
function, but we show how we can easily add alternate IDs given that we
know the organism and the source of the transcriptome. (This function
currently only works for Gencode and Ensembl gene or transcript IDs,
but could be extended to work for arbitrary sources.)


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
gse.int <- gse
assayNames(gse.int)
```

```
## [1] "counts"    "abundance" "length"
```

```r
# converting estimated counts to integer
assay(gse.int) <- round(assay(gse.int))
mode(assay(gse.int)) <- "integer"
# rename the "length" assay to "avgTxLength"
# DESeq2 will then use this matrix for calculating bias offsets
assayNames(gse.int)[3] <- "avgTxLength"
# turn the treatment variable into a factor
gse.int$treatment <- factor(gse.int$treatment)
# make a DESeqDataSet
dds <- DESeqDataSet(gse.int, ~treatment)
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
ah <- AnnotationHub()
```

```
## snapshotDate(): 2017-10-27
```

```r
chip <- query(ah, c("GM12878", "MEF2A", "narrowPeak"))[[1]]
```

```
## loading from cache '/Users/love//.AnnotationHub/28040'
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
comparing hg38 to hg38.


```r
url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
file <- "hg19ToHg38.over.chain.gz"
if (!file.exists(file)) download.file(url, file, method="wget")
system(paste0("gunzip ",file))
```

We move our ChIP-seq data to hg38:


```r
chip.lift <- liftOverHelper(chip, chainfile="hg19ToHg38.over.chain", to="hg38")
```

```
## Warning in var %in% extra[[i]]: closing unused connection 6
## (hg19ToHg38.over.chain)
```

Now we can find the nearest transcript to a given ChIP-seq peak:


```r
nearest(chip.lift[1], se)
```

```
## [1] 61708
```

```r
assay(se)[nearest(chip.lift[1], se),,drop=FALSE]
```

```
##                   SRR1565944 SRR1565945 SRR1565946 SRR1565947 SRR1565948
## ENST00000379350.5          0          0          0          0          0
##                   SRR1565949 SRR1565943 SRR1565939 SRR1565940 SRR1565942
## ENST00000379350.5          0          0          0          2          0
##                   SRR1565941 SRR1565938 SRR1565932 SRR1565933 SRR1565934
## ENST00000379350.5          1          0          0          0          0
##                   SRR1565935 SRR1565936 SRR1565937 SRR1565928 SRR1565926
## ENST00000379350.5          0          0          2          0          0
##                   SRR1565931 SRR1565930 SRR1565929 SRR1565927
## ENST00000379350.5          0          0          0          0
```

Or we can take a slice of the transcriptome data that is within 1
megabase of a given ChIP-seq peak:


```r
# which rows of SE in this window?
which(overlapsAny(se, chip.lift[1] + 1e6))
```

```
##  [1] 61702 61703 61704 61705 61706 61707 61708 61709 61710 61711 61712
## [12] 61713 61714 61715 61716 61717 61718 61719 61720 61721 61722 61723
## [23] 61724 61725 61726 61727 61728 61729 61730 61731 61732 61733 61734
## [34] 61735 61736 61737 61738 61739 61740 61741 61742 61743 61744 61745
## [45] 61746 61747 61748 61749 61750 61751 61752 61753 61754 61755 61756
## [56] 61757 61758 61759 61760 61761 61762 61763 61764 61765 61766 61767
## [67] 61768 61769 61770 61771 61772 61773 61774 61775 61776 61777 61778
## [78] 61779 61780 61781 61782 61783
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
##  $ salmon_version       : chr [1:24] "0.8.2" "0.8.2" "0.8.2" "0.8.2" ...
##  $ samp_type            : chr [1:24] "none" "none" "none" "none" ...
##  $ num_libraries        : int [1:24] 1 1 1 1 1 1 1 1 1 1 ...
##  $ library_types        : chr [1:24] "IU" "IU" "IU" "IU" ...
##  $ frag_dist_length     : int [1:24] 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 ...
##  $ seq_bias_correct     : logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ gc_bias_correct      : logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ num_bias_bins        : int [1:24] 4096 4096 4096 4096 4096 4096 4096 4096 4096 4096 ...
##  $ mapping_type         : chr [1:24] "mapping" "mapping" "mapping" "mapping" ...
##  $ num_targets          : int [1:24] 199324 199324 199324 199324 199324 199324 199324 199324 199324 199324 ...
##  $ serialized_eq_classes: logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ length_classes       : int [1:5, 1:24] 512 652 1007 2252 103053 512 652 1007 2252 103053 ...
##  $ index_seq_hash       : chr "13efe75909a197f8d30f7e17ed5f3f9f03a398b6796aa746416d74b4b6aa8aeb"
##  $ index_name_hash      : chr "59f80b32210db10dcd68ad4f8309b691b9c236659071dd50f9410e32e374f0a9"
##  $ num_bootstraps       : int [1:24] 0 0 0 0 0 0 0 0 0 0 ...
##  $ num_processed        : int [1:24] 6443960 6895927 6740421 6583297 6660199 6501138 6400538 6441883 6740421 6845343 ...
##  $ num_mapped           : int [1:24] 3352730 4172715 3433227 4108089 4428644 3711190 3331695 3411004 3446282 3887156 ...
##  $ percent_mapped       : num [1:24] 52 60.5 50.9 62.4 66.5 ...
##  $ call                 : chr [1:24] "quant" "quant" "quant" "quant" ...
##  $ start_time           : chr [1:24] "Wed Jun 21 10:45:14 2017" "Wed Jun 21 10:46:14 2017" "Wed Jun 21 10:47:20 2017" "Wed Jun 21 10:48:12 2017" ...
##  $ end_time             : chr [1:24] "Wed Jun 21 10:46:13 2017" "Wed Jun 21 10:47:20 2017" "Wed Jun 21 10:48:11 2017" "Wed Jun 21 10:49:16 2017" ...
```

```r
str(metadata(se)$txomeInfo)
```

```
## List of 7
##  $ index_seq_hash: chr "13efe75909a197f8d30f7e17ed5f3f9f03a398b6796aa746416d74b4b6aa8aeb"
##  $ source        : chr "Gencode"
##  $ organism      : chr "Homo sapiens"
##  $ version       : chr "26"
##  $ genome        : chr "GRCh38"
##  $ fasta         : chr "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz"
##  $ gtf           : chr "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz"
```

```r
str(metadata(se)$tximetaInfo)
```

```
## List of 2
##  $ version   :Classes 'package_version', 'numeric_version'  hidden list of 1
##   ..$ : int [1:3] 0 0 5
##  $ importTime: POSIXct[1:1], format: "2017-12-06 10:31:07"
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
##  version  R Under development (unstable) (2017-08-07 r73058)
##  system   x86_64, darwin15.6.0                              
##  ui       X11                                               
##  language (EN)                                              
##  collate  en_US.UTF-8                                       
##  tz       America/New_York                                  
##  date     2017-12-06
```

```
## Packages -----------------------------------------------------------------
```

```
##  package                           * version  date       source        
##  AnnotationDbi                     * 1.39.3   2017-09-11 Bioconductor  
##  AnnotationFilter                    1.1.9    2017-09-18 Bioconductor  
##  AnnotationHub                     * 2.9.19   2017-10-06 Bioconductor  
##  assertthat                          0.2.0    2017-04-11 CRAN (R 3.5.0)
##  backports                           1.1.1    2017-09-25 CRAN (R 3.5.0)
##  base                              * 3.5.0    2017-08-07 local         
##  bindr                               0.1      2016-11-13 CRAN (R 3.5.0)
##  bindrcpp                          * 0.2      2017-06-17 CRAN (R 3.5.0)
##  Biobase                           * 2.37.2   2017-08-08 Bioconductor  
##  BiocFileCache                     * 1.1.17   2017-10-06 Bioconductor  
##  BiocGenerics                      * 0.23.2   2017-10-06 Bioconductor  
##  BiocInstaller                     * 1.28.0   2017-11-07 Bioconductor  
##  BiocParallel                        1.11.9   2017-10-06 Bioconductor  
##  biomaRt                             2.33.4   2017-08-08 Bioconductor  
##  Biostrings                          2.45.4   2017-09-11 Bioconductor  
##  bit                                 1.1-12   2014-04-09 CRAN (R 3.5.0)
##  bit64                               0.9-7    2017-05-08 CRAN (R 3.5.0)
##  bitops                              1.0-6    2013-08-17 CRAN (R 3.5.0)
##  blob                                1.1.0    2017-06-17 CRAN (R 3.5.0)
##  BSgenome                            1.45.3   2017-09-18 Bioconductor  
##  codetools                           0.2-15   2016-10-05 CRAN (R 3.5.0)
##  commonmark                          1.4      2017-09-01 CRAN (R 3.5.0)
##  compiler                            3.5.0    2017-08-07 local         
##  crayon                              1.3.4    2017-09-16 CRAN (R 3.5.0)
##  curl                                3.0      2017-10-06 CRAN (R 3.5.0)
##  datasets                          * 3.5.0    2017-08-07 local         
##  DBI                                 0.7      2017-06-18 CRAN (R 3.5.0)
##  dbplyr                            * 1.1.0    2017-06-27 CRAN (R 3.5.0)
##  DelayedArray                      * 0.3.21   2017-10-06 Bioconductor  
##  desc                                1.1.1    2017-08-03 CRAN (R 3.5.0)
##  devtools                          * 1.13.3   2017-08-02 CRAN (R 3.5.0)
##  digest                              0.6.12   2017-01-27 CRAN (R 3.5.0)
##  dplyr                               0.7.4    2017-09-28 CRAN (R 3.5.0)
##  ensembldb                           2.1.13   2017-10-06 Bioconductor  
##  evaluate                            0.10.1   2017-06-24 CRAN (R 3.5.0)
##  GenomeInfoDb                      * 1.13.5   2017-10-06 Bioconductor  
##  GenomeInfoDbData                    0.99.1   2017-08-08 Bioconductor  
##  GenomicAlignments                   1.13.6   2017-09-18 Bioconductor  
##  GenomicFeatures                   * 1.29.11  2017-10-06 Bioconductor  
##  GenomicRanges                     * 1.29.15  2017-10-06 Bioconductor  
##  glue                                1.1.1    2017-06-21 CRAN (R 3.5.0)
##  GO.db                             * 3.4.1    2017-08-17 Bioconductor  
##  graph                               1.55.0   2017-08-22 Bioconductor  
##  graphics                          * 3.5.0    2017-08-07 local         
##  grDevices                         * 3.5.0    2017-08-07 local         
##  grid                                3.5.0    2017-08-07 local         
##  here                              * 0.1      2017-05-28 CRAN (R 3.5.0)
##  hms                                 0.3      2016-11-22 CRAN (R 3.5.0)
##  Homo.sapiens                      * 1.3.1    2017-12-06 Bioconductor  
##  htmltools                           0.3.6    2017-04-28 CRAN (R 3.5.0)
##  httpuv                              1.3.5    2017-07-04 CRAN (R 3.5.0)
##  httr                                1.3.1    2017-08-20 CRAN (R 3.5.0)
##  interactiveDisplayBase              1.15.0   2017-08-09 Bioconductor  
##  IRanges                           * 2.11.18  2017-10-06 Bioconductor  
##  jsonlite                            1.5      2017-06-01 CRAN (R 3.5.0)
##  knitr                               1.17     2017-08-10 CRAN (R 3.5.0)
##  lattice                             0.20-35  2017-03-25 CRAN (R 3.5.0)
##  lazyeval                            0.2.0    2016-06-12 CRAN (R 3.5.0)
##  magrittr                            1.5      2014-11-22 CRAN (R 3.5.0)
##  Matrix                              1.2-11   2017-08-16 CRAN (R 3.5.0)
##  matrixStats                       * 0.52.2   2017-04-14 CRAN (R 3.5.0)
##  memoise                             1.1.0    2017-04-21 CRAN (R 3.5.0)
##  methods                           * 3.5.0    2017-08-07 local         
##  mime                                0.5      2016-07-07 CRAN (R 3.5.0)
##  org.Hs.eg.db                      * 3.4.1    2017-09-07 Bioconductor  
##  OrganismDbi                       * 1.19.0   2017-09-07 Bioconductor  
##  parallel                          * 3.5.0    2017-08-07 local         
##  pkgconfig                           2.0.1    2017-03-21 CRAN (R 3.5.0)
##  prettyunits                         1.0.2    2015-07-13 CRAN (R 3.5.0)
##  progress                            1.1.2    2016-12-14 CRAN (R 3.5.0)
##  ProtGenerics                        1.9.1    2017-10-06 Bioconductor  
##  R6                                  2.2.2    2017-06-17 CRAN (R 3.5.0)
##  rappdirs                            0.3.1    2016-03-28 CRAN (R 3.5.0)
##  RBGL                                1.53.0   2017-08-22 Bioconductor  
##  Rcpp                                0.12.13  2017-09-28 CRAN (R 3.5.0)
##  RCurl                               1.95-4.8 2016-03-01 CRAN (R 3.5.0)
##  readr                             * 1.1.1    2017-05-16 CRAN (R 3.5.0)
##  rjson                               0.2.15   2014-11-03 CRAN (R 3.5.0)
##  rlang                               0.1.2    2017-08-09 CRAN (R 3.5.0)
##  rmarkdown                         * 1.6      2017-06-15 CRAN (R 3.5.0)
##  roxygen2                            6.0.1    2017-02-06 CRAN (R 3.5.0)
##  rprojroot                           1.2      2017-01-16 CRAN (R 3.5.0)
##  Rsamtools                           1.29.1   2017-09-11 Bioconductor  
##  RSQLite                             2.0      2017-06-19 CRAN (R 3.5.0)
##  rstudioapi                          0.7      2017-09-07 CRAN (R 3.5.0)
##  rtracklayer                       * 1.37.3   2017-08-08 Bioconductor  
##  S4Vectors                         * 0.15.11  2017-10-06 Bioconductor  
##  shiny                               1.0.5    2017-08-23 CRAN (R 3.5.0)
##  stats                             * 3.5.0    2017-08-07 local         
##  stats4                            * 3.5.0    2017-08-07 local         
##  stringi                             1.1.5    2017-04-07 CRAN (R 3.5.0)
##  stringr                             1.2.0    2017-02-18 CRAN (R 3.5.0)
##  SummarizedExperiment              * 1.7.10   2017-10-06 Bioconductor  
##  testthat                          * 1.0.2    2016-04-23 CRAN (R 3.5.0)
##  tibble                              1.3.4    2017-08-22 CRAN (R 3.5.0)
##  tools                               3.5.0    2017-08-07 local         
##  TxDb.Hsapiens.UCSC.hg19.knownGene * 3.2.2    2017-08-22 Bioconductor  
##  tximeta                           * 0.0.5    <NA>       Bioconductor  
##  tximport                            1.5.1    2017-10-06 Bioconductor  
##  utils                             * 3.5.0    2017-08-07 local         
##  withr                               2.0.0    2017-07-28 CRAN (R 3.5.0)
##  XML                                 3.98-1.9 2017-06-19 CRAN (R 3.5.0)
##  xml2                                1.1.1    2017-01-24 CRAN (R 3.5.0)
##  xtable                              1.8-2    2016-02-05 CRAN (R 3.5.0)
##  XVector                             0.17.1   2017-09-11 Bioconductor  
##  yaml                                2.1.14   2016-11-12 CRAN (R 3.5.0)
##  zlibbioc                            1.23.0   2017-08-08 Bioconductor
```
