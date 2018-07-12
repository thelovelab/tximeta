---
title: "tximeta: Working with derived transcriptomes"
author: "Michael Love, Rob Patro"
date: "07/12/2018"
output: 
  html_document:
    highlight: tango
---

# This package is in beta 

[See README](https://github.com/mikelove/tximeta/blob/master/README.md)



# Trying to import quants against a derived transcriptome


```r
dir <- system.file("extdata/derivedTxome/ERR188297", package="tximportData")
file <- file.path(dir, "quant.sf.gz")
file.exists(file)
```

```
## [1] TRUE
```

```r
coldata <- data.frame(files=file, names="SRR188297", sample="1",
                      stringsAsFactors=FALSE)
```

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
## couldn't find matching transcriptome, returning un-ranged SummarizedExperiment
```

# Documenting derived transcriptome for reproducible analysis

Clearly, for an arbitrary transcriptome used to generate the
quantification index, which does not directly match any
transcriptomes from known sources, there is not much that can be
done to populate the metadata during quantification import. However,
we can facilitate the following two cases: 

1) the transcriptome was derived locally from one or more publicly
available sources
2) the transcriptome was derived by another group, and they have
produced and shared a file that links their derived transcriptome
to publicly available sources

`tximeta` offers functionality to assist reproducible analysis in both
of these cases.

In the case of the quantification file above, the transcriptome was
generated locally by downloading the Ensembl cDNA FASTA for
Homo sapiens, version 87, and subsetting to the transcripts on the
"standard" chromsomes: 1-22, X, Y, MT. The following un-evaluated code
chunk reproduces the production of the transcriptome from publicly
available sources.


```r
library(Biostrings)
# Ensembl version 87
# ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
fa <- readDNAStringSet("Homo_sapiens.GRCh38.cdna.all.fa.gz")
chrs <- sapply(strsplit(sub(".* chromosome:GRCh38:(.*?) .*","\\1",names(fa)), ":"), `[`, 1)
chr.levels <- unique(sort(chrs))
idx <- chrs %in% c(1:22,"MT","X","Y")
sum(idx)
fa.sub <- fa[idx]
writeXStringSet(fa.sub, file="Homo_sapiens.GRCh38.cdna.std.chroms.fa")
```

To make this quantification reproducible, we make a
`derivedTxome` which records key information about the
source of the transcripts FASTA file, and the location of the
transcripts relative to a genome (GTF), it also records the signature
of the derived transcriptome that was computed by *Salmon* during the
`index` step. 

By default, `derivedTxome` will write out a JSON file
which can be shared with others, linking the signature of the index
with the other metadata, including FASTA and GTF sources. It will
write out to a file with the same name as the `indexDir`, but with
a `.json` extension added. This can be prevented with `write=FALSE`.


```r
dir <- system.file("extdata", package="tximeta")
indexDir <- file.path(dir, "Homo_sapiens.GRCh38.cdna.v87.std.chroms_salmon_0.8.1")
fastaFTP <- "ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
#gtfFTP <- "ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz"
gtfFTP <- file.path(dir,"Homo_sapiens.GRCh38.87.gtf.gz")
makeDerivedTxome(indexDir=indexDir,
                 source="Ensembl", organism="Homo sapiens",
                 version="87", genome="GRCh38",
                 fasta=fastaFTP, gtf=gtfFTP)
```

```
## writing derivedTxome to /Users/love/proj/tximeta/inst/extdata/Homo_sapiens.GRCh38.cdna.v87.std.chroms_salmon_0.8.1.json
```

```
## saving derivedTxome in bfc (first time)
```

After running `makeDerivedTxome`, the connection between this Salmon
index (and its signature) with the sources is saved for persistent
usage (this uses *BiocFileCache* and the specifics of the cache location
may change as `tximeta` develops).

With use of `tximeta` and a `derivedTxome` -- as with `tximeta` on a
standard, un-filtered transcriptome -- the software figures out if the
remote GTF has been accessed before, and on future calls, it will
simply load the stashed metadata (using *BiocFileCache*).


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
## found matching derived transcriptome:
## [ Ensembl - Homo sapiens - version 87 ]
## loading existing EnsDb created: 2018-07-12 12:49:51
## generating transcript ranges
```

We can see the appropriate metadata is attached (and here, using the
Ensembl-style notation of chromosomes, thanks to *ensembldb*).


```r
rowRanges(se)
```

```
## GRanges object with 162995 ranges and 6 metadata columns:
##                   seqnames              ranges strand |           tx_id
##                      <Rle>           <IRanges>  <Rle> |     <character>
##   ENST00000448914       14   22449113-22449125      + | ENST00000448914
##   ENST00000632684        7 142786213-142786224      + | ENST00000632684
##   ENST00000434970       14   22439007-22439015      + | ENST00000434970
##   ENST00000415118       14   22438547-22438554      + | ENST00000415118
##   ENST00000603693       15   21011451-21011469      - | ENST00000603693
##               ...      ...                 ...    ... .             ...
##   ENST00000454398        6   33131216-33143325      - | ENST00000454398
##   ENST00000457676       19   43303171-43304523      - | ENST00000457676
##   ENST00000457750        3 195614947-195620233      + | ENST00000457750
##   ENST00000520423        8   43539973-43540927      + | ENST00000520423
##   ENST00000523715        8   46549089-46550802      + | ENST00000523715
##                               tx_biotype tx_cds_seq_start tx_cds_seq_end
##                              <character>        <integer>      <integer>
##   ENST00000448914              TR_D_gene         22449113       22449125
##   ENST00000632684              TR_D_gene        142786213      142786224
##   ENST00000434970              TR_D_gene         22439007       22439015
##   ENST00000415118              TR_D_gene         22438547       22438554
##   ENST00000603693              IG_D_gene         21011451       21011469
##               ...                    ...              ...            ...
##   ENST00000454398 unprocessed_pseudogene             <NA>           <NA>
##   ENST00000457676 unprocessed_pseudogene             <NA>           <NA>
##   ENST00000457750 unprocessed_pseudogene             <NA>           <NA>
##   ENST00000520423 unprocessed_pseudogene             <NA>           <NA>
##   ENST00000523715 unprocessed_pseudogene             <NA>           <NA>
##                           gene_id         tx_name
##                       <character>     <character>
##   ENST00000448914 ENSG00000228985 ENST00000448914
##   ENST00000632684 ENSG00000282431 ENST00000632684
##   ENST00000434970 ENSG00000237235 ENST00000434970
##   ENST00000415118 ENSG00000223997 ENST00000415118
##   ENST00000603693 ENSG00000270451 ENST00000603693
##               ...             ...             ...
##   ENST00000454398 ENSG00000237398 ENST00000454398
##   ENST00000457676 ENSG00000230681 ENST00000457676
##   ENST00000457750 ENSG00000224769 ENST00000457750
##   ENST00000520423 ENSG00000253748 ENST00000520423
##   ENST00000523715 ENSG00000253425 ENST00000523715
##   -------
##   seqinfo: 47 sequences from GRCh38 genome
```

```r
seqinfo(se)
```

```
## Seqinfo object with 47 sequences from GRCh38 genome:
##   seqnames   seqlengths isCircular genome
##   1           248956422       <NA> GRCh38
##   10          133797422       <NA> GRCh38
##   11          135086622       <NA> GRCh38
##   12          133275309       <NA> GRCh38
##   13          114364328       <NA> GRCh38
##   ...               ...        ...    ...
##   KI270744.1     168472       <NA> GRCh38
##   KI270750.1     148850       <NA> GRCh38
##   MT              16569       <NA> GRCh38
##   X           156040895       <NA> GRCh38
##   Y            57227415       <NA> GRCh38
```

# Clear derivedTxomes (only for demonstration)


```r
library(BiocFileCache)
bfc <- BiocFileCache()
bfcinfo(bfc)
```

```
## # A tibble: 3 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
## 3 BFC9  deriv… 2018-07-12… 2018-07-12… /User… rela… 67f6…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

```r
bfcremove(bfc, bfcquery(bfc, "derivedTxomeDF")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 2 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

# Loading derivedTxome JSON files

If a collaborator or the Suppmentary Files for a publication 
shares a `derivedTxome` JSON file, we can likewise use `tximeta`
to automatically assemble the relevant metadata. This implies that the
other person has used `tximeta` with the function `makeDerivedTxome`
demonstrated above, pointing to their Salmon index and to the FASTA
and GTF source.

We simply point to the JSON file and use `loadDerivedTxome` and then
the relevant metadata is saved for persistent usage (using
*BiocFileCache*).


```r
dir <- system.file("extdata", package="tximeta")
jsonFile <- file.path(dir, "Homo_sapiens.GRCh38.cdna.v87.std.chroms_salmon_0.8.1.json")
loadDerivedTxome(jsonFile)
```

```
## saving derivedTxome in bfc (first time)
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
## found matching derived transcriptome:
## [ Ensembl - Homo sapiens - version 87 ]
## loading existing EnsDb created: 2018-07-12 12:49:51
## generating transcript ranges
```

# Clear derivedTxomes (only for demonstration)


```r
bfc <- BiocFileCache()
bfcinfo(bfc)
```

```
## # A tibble: 3 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
## 3 BFC10 deriv… 2018-07-12… 2018-07-12… /User… rela… 67f1…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

```r
bfcremove(bfc, bfcquery(bfc, "derivedTxomeDF")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 2 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC6  Homo_… 2018-07-12… 2018-07-12… /User… rela… 67f2…               NA
## 2 BFC7  genco… 2018-07-12… 2018-07-12… /User… rela… 67f7…               NA
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
##  date     2018-07-12
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version   date       source         
##  AnnotationDbi        * 1.42.1    2018-05-08 Bioconductor   
##  AnnotationFilter     * 1.4.0     2018-05-01 Bioconductor   
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
##  Biostrings           * 2.48.0    2018-05-01 Bioconductor   
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
##  desc                   1.2.0     2018-05-01 CRAN (R 3.5.0) 
##  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0) 
##  digest                 0.6.15    2018-01-28 cran (@0.6.15) 
##  dplyr                  0.7.5     2018-05-19 cran (@0.7.5)  
##  ensembldb            * 2.4.1     2018-05-07 Bioconductor   
##  evaluate               0.10.1    2017-06-24 CRAN (R 3.5.0) 
##  GenomeInfoDb         * 1.16.0    2018-05-01 Bioconductor   
##  GenomeInfoDbData       1.1.0     2018-01-10 Bioconductor   
##  GenomicAlignments      1.16.0    2018-05-01 Bioconductor   
##  GenomicFeatures      * 1.32.0    2018-05-01 Bioconductor   
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
##  rtracklayer            1.40.2    2018-05-08 Bioconductor   
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
##  tximeta              * 0.0.10    <NA>       Bioconductor   
##  tximport             * 1.8.0     2018-05-01 Bioconductor   
##  utf8                   1.1.4     2018-05-24 CRAN (R 3.5.0) 
##  utils                * 3.5.0     2018-04-24 local          
##  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0) 
##  XML                    3.98-1.11 2018-04-16 CRAN (R 3.5.0) 
##  xml2                   1.2.0     2018-01-24 CRAN (R 3.5.0) 
##  XVector              * 0.20.0    2018-05-01 Bioconductor   
##  yaml                   2.1.19    2018-05-01 CRAN (R 3.5.0) 
##  zlibbioc               1.26.0    2018-05-01 Bioconductor
```
