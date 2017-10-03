---
title: "tximeta: Import transcript quantification with automagic generation of metadata"
author: "Michael Love, Rob Patro"
date: "10/03/2017"
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
download.file("https://github.com/mikelove/asthma/archive/master.zip", dest)
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
```

```
## here() starts at /Users/love/proj/tximeta
```

```r
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
coldata$files <- here("extdata","asthma-master","data","quant",
                      coldata$run,"quant.sf.gz")
coldata$names <- coldata$run
```


```r
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(SummarizedExperiment))
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
## 1
```

```
## 2
```

```
## 3
```

```
## 4
```

```
## 5
```

```
## 6
```

```
## 7
```

```
## 8
```

```
## 9
```

```
## 10
```

```
## 11
```

```
## 12
```

```
## 13
```

```
## 14
```

```
## 15
```

```
## 16
```

```
## 17
```

```
## 18
```

```
## 19
```

```
## 20
```

```
## 21
```

```
## 22
```

```
## 23
```

```
## 24
```

```
## 
```

```
## found matching transcriptome:
## [ gencode - Homo sapiens - version 26 ]
```

```
## loading existing TxDb
```

```
## generating transcript ranges
```

```
## fetching genome info
```

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

Here we show the three matrices that were imported (note, this part
would need updating for un-reduced inferential variance matrices).


```r
assayNames(se)
```

```
## [1] "abundance" "counts"    "length"
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

# Demo: find nearest transcripts to a ChIP-seq peak

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
## snapshotDate(): 2017-08-31
```

```r
chip <- query(ah, c("GM12878", "MEF2A", "narrowPeak"))[[1]]
```

```
## require("rtracklayer")
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
if (!file.exists(file)) download.file(url, file)
system(paste0("gunzip ",file))
```

We move our ChIP-seq data to hg38:


```r
chip.lift <- liftOverHelper(chip, chainfile="hg19ToHg38.over.chain", to="hg38")
```

```
## Warning in .Internal(get(x, envir, mode, inherits)): closing unused
## connection 6 (hg19ToHg38.over.chain)
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
## ENST00000379350.5          0          0          0  0.0975893          0
##                   SRR1565941 SRR1565938 SRR1565932 SRR1565933 SRR1565934
## ENST00000379350.5  0.0493117          0          0          0          0
##                   SRR1565935 SRR1565936 SRR1565937 SRR1565928 SRR1565926
## ENST00000379350.5          0          0  0.0879775          0          0
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
## [1] "quantInfo"   "txomeInfo"   "tximetaInfo" "txdbInfo"
```

```r
str(metadata(se)$quantInfo)
```

```
## List of 21
##  $ salmon_version       : chr [1:24] "0.8.2" "0.8.2" "0.8.2" "0.8.2" ...
##  $ samp_type            : chr [1:24] "none" "none" "none" "none" ...
##  $ num_libraries        : num [1:24] 1 1 1 1 1 1 1 1 1 1 ...
##  $ library_types        : chr [1:24] "IU" "IU" "IU" "IU" ...
##  $ frag_dist_length     : num [1:24] 1001 1001 1001 1001 1001 ...
##  $ seq_bias_correct     : logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ gc_bias_correct      : logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ num_bias_bins        : num [1:24] 4096 4096 4096 4096 4096 ...
##  $ mapping_type         : chr [1:24] "mapping" "mapping" "mapping" "mapping" ...
##  $ num_targets          : num [1:24] 199324 199324 199324 199324 199324 ...
##  $ serialized_eq_classes: logi [1:24] FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ length_classes       : num [1:5, 1:24] 512 652 1007 2252 103053 ...
##  $ index_seq_hash       : chr "13efe75909a197f8d30f7e17ed5f3f9f03a398b6796aa746416d74b4b6aa8aeb"
##  $ index_name_hash      : chr "59f80b32210db10dcd68ad4f8309b691b9c236659071dd50f9410e32e374f0a9"
##  $ num_bootstraps       : num [1:24] 0 0 0 0 0 0 0 0 0 0 ...
##  $ num_processed        : num [1:24] 6443960 6895927 6740421 6583297 6660199 ...
##  $ num_mapped           : num [1:24] 3352730 4172715 3433227 4108089 4428644 ...
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
##  $ source        : chr "gencode"
##  $ organism      : chr "Homo sapiens"
##  $ version       : int 26
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
##   ..$ : int [1:3] 0 0 1
##  $ importTime: POSIXct[1:1], format: "2017-09-01 15:37:05"
```

```r
str(metadata(se)$txdbInfo)
```

```
##  Named chr [1:15] "TxDb" "GenomicFeatures" ...
##  - attr(*, "names")= chr [1:15] "Db type" "Supporting package" "Data source" "Organism" ...
```

# Next steps

### Challenges

* Building out actual, sustainable plan for supporting as many
  organisms and sources as possible. We can define rules which
  determine where the FASTA and GTF files will be based on `source` and
  `version` (also here I ignored something like "type", e.g. CHR
  or ALL gene files from Gencode)
* Want to support computational reproducibility for "derived
  transcriptomes" or de novo transcriptomes. See this 
  [GitHub Issue](https://github.com/mikelove/tximeta/issues/2).
  How can we encapsulate the steps to reproduce the generation of
  the transcriptome from its public source, and use signature to
  guarantee we have the same sequence.
* Facilitate functional annotation, either with vignettes/workflow or
  with additional functionality. E.g.: 
  housekeeping genes, arbitrary gene sets, genes expressed in GTEx tissues
* liftOver is clunky and doesn't integrate with
  GenomeInfoDb. It requires user input and there's a chance to
  mis-annotate. Ideally this should all be automated.

### Facilitate plots and summaries
    
* Basic plots across samples: abundances, mapping rates, rich bias model parameters
* Time summaries: when quantified? when imported? I would love to
  know when the library was prepared and sequenced but this seems hopeless.

### Basic functionality

* Switching `rowRanges` from transcript ranges to exons-by-transcript ranges list
* Summarization to gene-level (where building `tx2gene` is no longer necessary)
* As is already supported in the `tximport` release, also import inferential
  variance matrices (Gibbs samples or bootstrap samples)

