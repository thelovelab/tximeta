---
title: "tximeta: Working with linked transcriptomes"
author: "Michael Love, Rob Patro, Charlotte Soneson, Peter Hickey"
date: "07/16/2018"
output: 
  rmarkdown::html_document:
    highlight: tango
abstract: >
  `tximeta` automatically imports relevant metadata when the
  transcriptome matches a known source, but also facilitates the
  linking of transcriptomes used as for a *Salmon* index with relevant
  public sources. The linking is important in the case that the
  transcript sequence no longer matches a known source (combined or
  filtered FASTA files), or if the source is not known to
  `tximeta`. Here we demonstrate how to make a *linkedTxome* and how
  to share and load a *linkedTxome*.
vignette: |
  %\VignetteIndexEntry{Working with linked transcriptomes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Trying to import quantification files with an unknown transcriptome

Here we point to *Salmon* quantification files which were quantified
against a transcriptome combining two Ensembl FASTA files: the cDNA
and the non-coding transcripts for *Drosophila melanogaster*.


```r
dir <- system.file("extdata/salmon_dm/SRR1197474", package="tximportData")
file <- file.path(dir, "quant.sf.gz")
file.exists(file)
```

```
## [1] TRUE
```

```r
coldata <- data.frame(files=file, names="SRR1197474", sample="1",
                      stringsAsFactors=FALSE)
```

Trying to import the files gives a message that `tximeta` couldn't find
a matching transcriptome, so it returns an un-ranged
*SummarizedExperiment*. 


```r
suppressPackageStartupMessages(library(tximeta))
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
## couldn't find matching transcriptome, returning un-ranged SummarizedExperiment
```

# Documenting linked transcriptome for reproducible analysis

If the transcriptome used to generate the *Salmon* index does not
match any transcriptomes from known sources (e.g. from combining or filtering
known transcriptome files), there is not much that can be done to
automatically populate the metadata during quantification
import. However, we can facilitate the following two cases: 

1) the transcriptome was created locally and has been linked to its
public source(s) 
2) the transcriptome was produced by another group, and
they have produced and shared a file that links the transcriptome to
public source(s)

`tximeta` offers functionality to assist reproducible analysis in both
of these cases.

In the case of the quantification file above, the transcriptome was
generated locally by downloading and combining the Ensembl cDNA and
non-coding FASTA files *Drosophila melanogaster*, version 92. The
following un-evaluated command line code chunk reproduces the
production of the transcriptome from publicly available sources.

```
wget ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.ncrna.fa.gz
cat Drosophila_melanogaster.BDGP6.cdna.all.fa.gz Drosophila_melanogaster.BDGP6.ncrna.fa.gz > Drosophila_melanogaster.BDGP6.v92.fa.gz
```

To make this quantification reproducible, we make a `linkedTxome`
which records key information about the sources of the transcript
FASTA files, and the location of the relevant GTF file. It also
records the signature of the transcriptome that was computed by
*Salmon* during the `index` step.

By default, `linkedTxome` will write out a JSON file which can be
shared with others, linking the signature of the index with the other
metadata, including FASTA and GTF sources. By default, it will write
out to a file with the same name as the `indexDir`, but with a `.json`
extension added. This can be prevented with `write=FALSE`, and the
file location can be changed with `jsonFile`.

First we specify the path where the *Salmon* index is located. 

**Note**: typically you would not use `system.file` to find this
directory, but simply define `indexDir` to be the path of the
*Salmon* directory on your machine. Here we use `system.file` because
we have included parts of a *Salmon* index directory in the *tximeta*
package itself for demonstration of functionality in this vignette.


```r
dir <- system.file("extdata", package="tximeta")
indexDir <- file.path(dir, "Drosophila_melanogaster.BDGP6.v92_salmon_0.10.2")
```

Now we provide the location of the FASTA files and the GTF file for
this transcriptome.


```r
fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz",
              "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.ncrna.fa.gz")
#gtfFTP <- "ftp://ftp.ensembl.org/pub/release-92/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.92.gtf.gz"
```

Instead of the above commented-out FTP location for the GTF file, we
specify a location within an R package. This step is just to avoid
downloading from a remote FTP during vignette building. This use of
`system.file` to point to a file in an R package is specific to this
vignette and would not be used in a typical workflow.


```r
dir2 <- system.file("extdata/salmon_dm", package="tximportData")
gtfPath <- file.path(dir2,"Drosophila_melanogaster.BDGP6.92.gtf.gz")
```

Finally, we create a *linkedTxome*.  In this vignette, we point to a
temporary directory for the JSON file, but a more typical workflow
would write the JSON file to the same location as the *Salmon* index
by not specifying `jsonFile`.

`makeLinkedTxome` performs two operation: (1) it creates a new entry in
an internal table that links the transcriptome used in the *Salmon*
index to its sources, and (2) it creates a JSON file such that this
*linkedTxome* can be shared.


```r
tmp <- tempdir()
jsonFile <- file.path(tmp, paste0(basename(indexDir), ".json"))
makeLinkedTxome(indexDir=indexDir,
                 source="Ensembl", organism="Drosophila melanogaster",
                 version="92", genome="BDGP6",
                 fasta=fastaFTP, gtf=gtfPath,
                 jsonFile=jsonFile)
```

```
## writing linkedTxome to /var/folders/kd/xns8y9c51sz668rrjn0bvlqh0000gn/T//Rtmpj6yQ9O/Drosophila_melanogaster.BDGP6.v92_salmon_0.10.2.json
```

```
## saving linkedTxome in bfc (first time)
```

After running `makeLinkedTxome`, the connection between this *Salmon*
index (and its signature) with the sources is saved for persistent
usage.

With use of `tximeta` and a `linkedTxome` -- as with `tximeta` on a
standard, un-filtered transcriptome -- the software figures out if the
remote GTF has been accessed before, and on future calls, it will
simply load the pre-computed metadata and transcript ranges.

Note the warning that 9 of the transcripts are missing from the GTF
file and so are dropped from the final output. This is a problem
coming from the annotation source, and not easily avoided by
`tximeta`. 


```r
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
##  9 out of 33681 are missing from the GTF and dropped from SummarizedExperiment output
```

We can see that the appropriate metadata and transcript ranges are
attached.


```r
rowRanges(se)
```

```
## GRanges object with 33672 ranges and 6 metadata columns:
##               seqnames            ranges strand |       tx_id
##                  <Rle>         <IRanges>  <Rle> | <character>
##   FBtr0075502       3L 15808322-15808883      + | FBtr0075502
##   FBtr0300738       3R   5783105-5787336      + | FBtr0300738
##   FBtr0300739       3R   5781762-5787336      + | FBtr0300739
##   FBtr0300737       3R   5781762-5787336      + | FBtr0300737
##   FBtr0300736       3R   5783105-5787336      + | FBtr0300736
##           ...      ...               ...    ... .         ...
##   FBtr0086850       2R 17701229-17701297      + | FBtr0086850
##   FBtr0113576       3R   5596201-5596340      - | FBtr0113576
##   FBtr0076635       3L   8601948-8602031      + | FBtr0076635
##   FBtr0309760       3L     891250-891475      + | FBtr0309760
##   FBtr0113549       2L 20419932-20420065      + | FBtr0113549
##                   tx_biotype tx_cds_seq_start tx_cds_seq_end     gene_id
##                  <character>        <integer>      <integer> <character>
##   FBtr0075502 protein_coding         15808418       15808716 FBgn0036531
##   FBtr0300738 protein_coding          5783217        5787117 FBgn0037375
##   FBtr0300739 protein_coding          5781900        5787117 FBgn0037375
##   FBtr0300737 protein_coding          5781900        5787117 FBgn0037375
##   FBtr0300736 protein_coding          5783217        5787117 FBgn0037375
##           ...            ...              ...            ...         ...
##   FBtr0086850         snoRNA             <NA>           <NA> FBgn0063388
##   FBtr0113576         snoRNA             <NA>           <NA> FBgn0082961
##   FBtr0076635         snoRNA             <NA>           <NA> FBgn0060291
##   FBtr0309760         snoRNA             <NA>           <NA> FBgn0263461
##   FBtr0113549         snoRNA             <NA>           <NA> FBgn0083032
##                   tx_name
##               <character>
##   FBtr0075502 FBtr0075502
##   FBtr0300738 FBtr0300738
##   FBtr0300739 FBtr0300739
##   FBtr0300737 FBtr0300737
##   FBtr0300736 FBtr0300736
##           ...         ...
##   FBtr0086850 FBtr0086850
##   FBtr0113576 FBtr0113576
##   FBtr0076635 FBtr0076635
##   FBtr0309760 FBtr0309760
##   FBtr0113549 FBtr0113549
##   -------
##   seqinfo: 25 sequences from BDGP6 genome
```

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

# Clear *linkedTxomes*

The following code removes the table with information about the
*linkedTxomes*. This is just for demonstration, so that we can show
how to load a JSON file below. 

**Note:** Running this code will clear information about any
*linkedTxomes* on your machine.


```r
library(BiocFileCache)
bfcloc <- getTximetaBFC()
bfc <- BiocFileCache(bfcloc)
bfcinfo(bfc)
```

```
## # A tibble: 4 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
## 3 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## 4 BFC33 linke… 2018-07-14… 2018-07-14… /User… rela… 36b5…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

```r
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 3 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
## 3 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

# Loading *linkedTxome* JSON files

If a collaborator or the Suppmentary Files for a publication shares a
`linkedTxome` JSON file, we can likewise use `tximeta` to
automatically assemble the relevant metadata and transcript
ranges. This implies that the other person has used `tximeta` with the
function `makeLinkedTxome` demonstrated above, pointing to their
*Salmon* index and to the FASTA and GTF source.

We point to the JSON file and use `loadLinkedTxome` and then the
relevant metadata is saved for persistent usage. In this case, we
saved the JSON file in a temporary directory.


```r
jsonFile <- file.path(tmp, paste0(basename(indexDir), ".json"))
loadLinkedTxome(jsonFile)
```

```
## saving linkedTxome in bfc (first time)
```

Again, using `tximeta` figures out whether it needs to access the
remote GTF or not, and assembles the appropriate object on the user's
behalf.


```r
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
##  9 out of 33681 are missing from the GTF and dropped from SummarizedExperiment output
```

# Clear *linkedTxomes*


```r
bfcloc <- getTximetaBFC()
bfc <- BiocFileCache(bfcloc)
bfcinfo(bfc)
```

```
## # A tibble: 4 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
## 3 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## 4 BFC34 linke… 2018-07-14… 2018-07-14… /User… rela… 36b5…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

```r
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 3 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
## 3 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

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
##  package              * version   date       source         
##  AnnotationDbi          1.42.1    2018-05-08 Bioconductor   
##  AnnotationFilter       1.4.0     2018-05-01 Bioconductor   
##  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0) 
##  backports              1.1.2     2017-12-13 cran (@1.1.2)  
##  base                 * 3.5.0     2018-04-24 local          
##  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0) 
##  bindrcpp             * 0.2.2     2018-03-29 CRAN (R 3.5.0) 
##  Biobase              * 2.40.0    2018-05-01 Bioconductor   
##  BiocFileCache        * 1.4.0     2018-05-01 Bioconductor   
##  BiocGenerics         * 0.26.0    2018-05-01 Bioconductor   
##  BiocInstaller        * 1.30.0    2018-05-04 Bioconductor   
##  BiocParallel         * 1.14.1    2018-05-06 Bioconductor   
##  biomaRt                2.36.1    2018-05-24 Bioconductor   
##  Biostrings             2.48.0    2018-05-01 Bioconductor   
##  bit                    1.1-14    2018-05-29 CRAN (R 3.5.0) 
##  bit64                  0.9-7     2017-05-08 CRAN (R 3.5.0) 
##  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0) 
##  blob                   1.1.1     2018-03-25 CRAN (R 3.5.0) 
##  cli                    1.0.0     2017-11-05 CRAN (R 3.5.0) 
##  codetools              0.2-15    2016-10-05 CRAN (R 3.5.0) 
##  commonmark             1.5       2018-04-28 CRAN (R 3.5.0) 
##  compiler               3.5.0     2018-04-24 local          
##  crayon                 1.3.4     2017-09-16 CRAN (R 3.5.0) 
##  curl                   3.2       2018-03-28 CRAN (R 3.5.0) 
##  datasets             * 3.5.0     2018-04-24 local          
##  DBI                    1.0.0     2018-05-02 CRAN (R 3.5.0) 
##  dbplyr               * 1.2.1     2018-02-19 CRAN (R 3.5.0) 
##  DelayedArray         * 0.6.1     2018-06-15 Bioconductor   
##  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0) 
##  digest                 0.6.15    2018-01-28 cran (@0.6.15) 
##  dplyr                  0.7.5     2018-05-19 cran (@0.7.5)  
##  ensembldb              2.4.1     2018-05-07 Bioconductor   
##  evaluate               0.10.1    2017-06-24 CRAN (R 3.5.0) 
##  GenomeInfoDb         * 1.16.0    2018-05-01 Bioconductor   
##  GenomeInfoDbData       1.1.0     2018-01-10 Bioconductor   
##  GenomicAlignments      1.16.0    2018-05-01 Bioconductor   
##  GenomicFeatures        1.32.0    2018-05-01 Bioconductor   
##  GenomicRanges        * 1.32.3    2018-05-16 Bioconductor   
##  glue                   1.2.0     2017-10-29 CRAN (R 3.5.0) 
##  graphics             * 3.5.0     2018-04-24 local          
##  grDevices            * 3.5.0     2018-04-24 local          
##  grid                   3.5.0     2018-04-24 local          
##  hms                    0.4.2     2018-03-10 CRAN (R 3.5.0) 
##  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0) 
##  httr                   1.3.1     2017-08-20 CRAN (R 3.5.0) 
##  IRanges              * 2.14.10   2018-05-16 Bioconductor   
##  jsonlite               1.5       2017-06-01 CRAN (R 3.5.0) 
##  knitr                  1.20      2018-02-20 CRAN (R 3.5.0) 
##  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0) 
##  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0) 
##  magrittr               1.5       2014-11-22 CRAN (R 3.5.0) 
##  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0) 
##  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.5.0) 
##  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0) 
##  methods              * 3.5.0     2018-04-24 local          
##  parallel             * 3.5.0     2018-04-24 local          
##  pillar                 1.2.3     2018-05-25 CRAN (R 3.5.0) 
##  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0) 
##  prettyunits            1.0.2     2015-07-13 CRAN (R 3.5.0) 
##  progress               1.2.0     2018-06-14 CRAN (R 3.5.0) 
##  ProtGenerics           1.12.0    2018-05-01 Bioconductor   
##  purrr                  0.2.5     2018-05-29 cran (@0.2.5)  
##  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0) 
##  rappdirs               0.3.1     2016-03-28 CRAN (R 3.5.0) 
##  Rcpp                   0.12.17   2018-05-18 cran (@0.12.17)
##  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.5.0) 
##  readr                  1.1.1     2017-05-16 CRAN (R 3.5.0) 
##  rlang                  0.2.1     2018-05-30 cran (@0.2.1)  
##  rmarkdown            * 1.9       2018-03-01 CRAN (R 3.5.0) 
##  roxygen2               6.0.1     2017-02-06 CRAN (R 3.5.0) 
##  rprojroot              1.3-2     2018-01-03 cran (@1.3-2)  
##  Rsamtools              1.32.2    2018-07-03 Bioconductor   
##  RSQLite                2.1.1     2018-05-06 CRAN (R 3.5.0) 
##  rtracklayer            1.40.3    2018-06-02 cran (@1.40.3) 
##  S4Vectors            * 0.18.3    2018-06-08 Bioconductor   
##  stats                * 3.5.0     2018-04-24 local          
##  stats4               * 3.5.0     2018-04-24 local          
##  stringi                1.2.3     2018-06-12 CRAN (R 3.5.0) 
##  stringr                1.3.1     2018-05-10 CRAN (R 3.5.0) 
##  SummarizedExperiment * 1.10.1    2018-05-11 Bioconductor   
##  testthat             * 2.0.0     2017-12-13 CRAN (R 3.5.0) 
##  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0) 
##  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0) 
##  tools                  3.5.0     2018-04-24 local          
##  tximeta              * 0.0.14    <NA>       Bioconductor   
##  tximport               1.8.0     2018-05-01 Bioconductor   
##  utf8                   1.1.4     2018-05-24 CRAN (R 3.5.0) 
##  utils                * 3.5.0     2018-04-24 local          
##  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0) 
##  XML                    3.98-1.11 2018-04-16 CRAN (R 3.5.0) 
##  xml2                   1.2.0     2018-01-24 CRAN (R 3.5.0) 
##  XVector                0.20.0    2018-05-01 Bioconductor   
##  yaml                   2.1.19    2018-05-01 CRAN (R 3.5.0) 
##  zlibbioc               1.26.0    2018-05-01 Bioconductor
```
