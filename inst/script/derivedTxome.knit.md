---
title: "tximeta: Working with derived transcriptomes"
author: "Michael Love, Rob Patro"
date: "12/04/2017"
output: 
  html_document:
    highlight: tango
---

# This package is in beta 

[See README](https://github.com/mikelove/tximeta/blob/master/README.md)

# Setup

First, to try out `tximeta` you'll need the example data, which is
contained in this GitHub repo. Here we download the whole repo as a
ZIP file.




```r
library(here)
dest <- here("extdata","crohns.zip")
download.file("https://github.com/mikelove/crohns/archive/master.zip", dest, method="wget")
unzip(dest, exdir=here("extdata"))
```

# Trying to import quants against a derived transcriptome


```r
library(here)
file <- here("extdata/crohns-master/data/ensembl/SRR1813877/quant.sf.gz")
coldata <- data.frame(sample="1", files=file, names="SRR1813877",
                      stringsAsFactors=FALSE)
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
produced and shared an object that links their derived transcriptome
to publicly available sources

`tximeta` offers functionality to assist reproducible analysis in both
of these cases.

In the case of the quantification file above, the transcriptome was
generated locally by downloading the Ensembl cDNA FASTA and GTF for
Homo sapiens, version 90, and subsetting to the transcripts on the
"standard" chromsomes: 1-22, X, Y, MT. The following code reproduces
the production of the transcriptome from publicly available sources:


```r
library(ensembldb)
dbfile <- ensDbFromGtf("Homo_sapiens.GRCh38.90.gtf.gz")
edb <- EnsDb("Homo_sapiens.GRCh38.90.sqlite")
# this is similar to Gencode CHR (except gencode has ~150 additional transcripts in the PAR of Y)
txps <- transcripts(edb, filter=AnnotationFilterList(SeqNameFilter(c(1:22, "X", "Y","MT"))))
library(Biostrings)
# this is missing some of the transcripts in the above
seqs <- readDNAStringSet("Homo_sapiens.GRCh38.cdna.all.fa.gz")
names(seqs) <- sub("(.*?)\\..*","\\1",names(seqs)) # split at '.'
common <- intersect(names(txps), names(seqs)) 
seqs.sub <- seqs[common]
writeXStringSet(seqs.sub, filepath="Homo_sapiens.GRCh38.cdna.std.chroms.fa")
```

To make this quantification reproducible, we make a
`derivedTxome` object which records key information about the
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
indexDir <- "Homo_sapiens.GRCh38.cdna.std.chroms_salmon_0.8.2"
fastaFTP <- "ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
gtfFTP <- "ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz"
makeDerivedTxome(indexDir=indexDir,
                 source="Ensembl", organism="Homo sapiens",
                 version="90", genome="GRCh38",
                 fasta=fastaFTP, gtf=gtfFTP)
```

```
## writing derivedTxome to Homo_sapiens.GRCh38.cdna.std.chroms_salmon_0.8.2.json
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
## [ Ensembl - Homo sapiens - version 90 ]
## loading existing EnsDb created: 2017-12-04 20:17:48
## generating transcript ranges
```

We can see the appropriate metadata is attached (and here, using the
Ensembl-style notation of chromosomes, thanks to *ensembldb*).


```r
rowRanges(se)
```

```
## GRanges object with 164776 ranges and 6 metadata columns:
##                   seqnames               ranges strand |           tx_id
##                      <Rle>            <IRanges>  <Rle> |     <character>
##   ENST00000456328        1       [11869, 14409]      + | ENST00000456328
##   ENST00000450305        1       [12010, 13670]      + | ENST00000450305
##   ENST00000488147        1       [14404, 29570]      - | ENST00000488147
##   ENST00000606857        1       [52473, 53312]      + | ENST00000606857
##   ENST00000642116        1       [57598, 64116]      + | ENST00000642116
##               ...      ...                  ...    ... .             ...
##   ENST00000420810        Y [26549425, 26549743]      + | ENST00000420810
##   ENST00000456738        Y [26586642, 26591601]      - | ENST00000456738
##   ENST00000435945        Y [26594851, 26634652]      - | ENST00000435945
##   ENST00000435741        Y [26626520, 26627159]      - | ENST00000435741
##   ENST00000431853        Y [56855244, 56855488]      + | ENST00000431853
##                                           tx_biotype tx_cds_seq_start
##                                          <character>        <integer>
##   ENST00000456328               processed_transcript             <NA>
##   ENST00000450305 transcribed_unprocessed_pseudogene             <NA>
##   ENST00000488147             unprocessed_pseudogene             <NA>
##   ENST00000606857             unprocessed_pseudogene             <NA>
##   ENST00000642116               processed_transcript             <NA>
##               ...                                ...              ...
##   ENST00000420810               processed_pseudogene             <NA>
##   ENST00000456738             unprocessed_pseudogene             <NA>
##   ENST00000435945             unprocessed_pseudogene             <NA>
##   ENST00000435741               processed_pseudogene             <NA>
##   ENST00000431853               processed_pseudogene             <NA>
##                   tx_cds_seq_end         gene_id         tx_name
##                        <integer>     <character>     <character>
##   ENST00000456328           <NA> ENSG00000223972 ENST00000456328
##   ENST00000450305           <NA> ENSG00000223972 ENST00000450305
##   ENST00000488147           <NA> ENSG00000227232 ENST00000488147
##   ENST00000606857           <NA> ENSG00000268020 ENST00000606857
##   ENST00000642116           <NA> ENSG00000240361 ENST00000642116
##               ...            ...             ...             ...
##   ENST00000420810           <NA> ENSG00000224240 ENST00000420810
##   ENST00000456738           <NA> ENSG00000227629 ENST00000456738
##   ENST00000435945           <NA> ENSG00000237917 ENST00000435945
##   ENST00000435741           <NA> ENSG00000231514 ENST00000435741
##   ENST00000431853           <NA> ENSG00000235857 ENST00000431853
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
bfc <- BiocFileCache(".")
bfcinfo(bfc)
```

```
## # A tibble: 3 x 8
##     rid                         rname         create_time
##   <chr>                         <chr>               <chr>
## 1 BFC20 Homo_sapiens.GRCh38.90.gtf.gz 2017-12-04 20:17:48
## 2 BFC22 gencode.v26.annotation.gtf.gz 2017-12-04 21:04:56
## 3 BFC23                derivedTxomeDF 2017-12-04 21:41:22
## # ... with 5 more variables: access_time <chr>, rpath <chr>, rtype <chr>,
## #   fpath <chr>, last_modified_time <chr>
```

```r
bfcremove(bfc, bfcquery(bfc, "derivedTxomeDF")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 2 x 8
##     rid                         rname         create_time
##   <chr>                         <chr>               <chr>
## 1 BFC20 Homo_sapiens.GRCh38.90.gtf.gz 2017-12-04 20:17:48
## 2 BFC22 gencode.v26.annotation.gtf.gz 2017-12-04 21:04:56
## # ... with 5 more variables: access_time <chr>, rpath <chr>, rtype <chr>,
## #   fpath <chr>, last_modified_time <chr>
```

# Loading derivedTxome JSON files

If a collaborator or the Suppmentary Files for a publication 
shares a `derivedTxomes` JSON file, we can likewise use `tximeta`
to automatically assemble the relevant metadata. This implies that the
other person has used `tximeta` with the function `makeDerivedTxome`
demonstrated above, pointing to their Salmon index and to the FASTA
and GTF source.

We simply point to the JSON file and use `loadDerivedTxome` and then
the relevant metadata is saved for persistent usage (using
*BiocFileCache*).


```r
jsonfile <- "Homo_sapiens.GRCh38.cdna.std.chroms_salmon_0.8.2.json"
loadDerivedTxome(jsonfile)
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
## [ Ensembl - Homo sapiens - version 90 ]
## loading existing EnsDb created: 2017-12-04 20:17:48
## generating transcript ranges
```

# Clear derivedTxomes (only for demonstration)


```r
bfc <- BiocFileCache(".")
bfcinfo(bfc)
```

```
## # A tibble: 3 x 8
##     rid                         rname         create_time
##   <chr>                         <chr>               <chr>
## 1 BFC20 Homo_sapiens.GRCh38.90.gtf.gz 2017-12-04 20:17:48
## 2 BFC22 gencode.v26.annotation.gtf.gz 2017-12-04 21:04:56
## 3 BFC24                derivedTxomeDF 2017-12-04 21:41:58
## # ... with 5 more variables: access_time <chr>, rpath <chr>, rtype <chr>,
## #   fpath <chr>, last_modified_time <chr>
```

```r
bfcremove(bfc, bfcquery(bfc, "derivedTxomeDF")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 2 x 8
##     rid                         rname         create_time
##   <chr>                         <chr>               <chr>
## 1 BFC20 Homo_sapiens.GRCh38.90.gtf.gz 2017-12-04 20:17:48
## 2 BFC22 gencode.v26.annotation.gtf.gz 2017-12-04 21:04:56
## # ... with 5 more variables: access_time <chr>, rpath <chr>, rtype <chr>,
## #   fpath <chr>, last_modified_time <chr>
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
##  version  R Under development (unstable) (2017-05-23 r72721)
##  system   x86_64, linux-gnu                                 
##  ui       X11                                               
##  language en_US                                             
##  collate  en_US.UTF-8                                       
##  tz       posixrules                                        
##  date     2017-12-04
```

```
## Packages -----------------------------------------------------------------
```

```
##  package                     * version  date       source        
##  AnnotationDbi               * 1.39.4   2017-10-18 Bioconductor  
##  AnnotationFilter              1.1.9    2017-10-08 Bioconductor  
##  AnnotationHub               * 2.9.19   2017-10-08 Bioconductor  
##  assertthat                    0.2.0    2017-04-11 CRAN (R 3.5.0)
##  backports                     1.1.1    2017-09-25 CRAN (R 3.5.0)
##  base                        * 3.5.0    2017-05-24 local         
##  bindr                         0.1      2016-11-13 CRAN (R 3.5.0)
##  bindrcpp                    * 0.2      2017-06-17 CRAN (R 3.5.0)
##  Biobase                     * 2.37.2   2017-05-24 Bioconductor  
##  BiocFileCache               * 1.2.0    2017-12-01 Bioconductor  
##  BiocGenerics                * 0.23.3   2017-10-08 Bioconductor  
##  BiocInstaller               * 1.28.0   2017-11-08 Bioconductor  
##  BiocParallel                  1.11.11  2017-10-18 Bioconductor  
##  biomaRt                       2.33.4   2017-10-08 Bioconductor  
##  Biostrings                  * 2.45.4   2017-10-08 Bioconductor  
##  bit                           1.1-12   2014-04-09 CRAN (R 3.5.0)
##  bit64                         0.9-7    2017-05-08 CRAN (R 3.5.0)
##  bitops                        1.0-6    2013-08-17 CRAN (R 3.5.0)
##  blob                          1.1.0    2017-06-17 CRAN (R 3.5.0)
##  BSgenome                    * 1.45.3   2017-10-08 Bioconductor  
##  BSgenome.Hsapiens.UCSC.hg19 * 1.4.0    2017-09-08 Bioconductor  
##  codetools                     0.2-15   2016-10-05 CRAN (R 3.5.0)
##  commonmark                    1.4      2017-09-01 CRAN (R 3.5.0)
##  compiler                      3.5.0    2017-05-24 local         
##  crayon                        1.3.4    2017-09-16 CRAN (R 3.5.0)
##  curl                          3.0      2017-10-06 CRAN (R 3.5.0)
##  datasets                    * 3.5.0    2017-05-24 local         
##  DBI                           0.7      2017-06-18 CRAN (R 3.5.0)
##  dbplyr                      * 1.1.0    2017-06-27 CRAN (R 3.5.0)
##  DelayedArray                * 0.3.21   2017-10-08 Bioconductor  
##  desc                          1.1.1    2017-08-03 CRAN (R 3.5.0)
##  devtools                    * 1.13.3   2017-08-02 CRAN (R 3.5.0)
##  digest                        0.6.12   2017-01-27 CRAN (R 3.5.0)
##  dplyr                         0.7.4    2017-09-28 CRAN (R 3.5.0)
##  ensembldb                     2.1.14   2017-10-18 Bioconductor  
##  evaluate                      0.10.1   2017-06-24 CRAN (R 3.5.0)
##  GenomeInfoDb                * 1.13.5   2017-10-08 Bioconductor  
##  GenomeInfoDbData              0.99.1   2017-10-08 Bioconductor  
##  GenomicAlignments             1.13.6   2017-10-08 Bioconductor  
##  GenomicFeatures             * 1.29.13  2017-10-18 Bioconductor  
##  GenomicRanges               * 1.29.15  2017-10-08 Bioconductor  
##  glue                          1.1.1    2017-06-21 CRAN (R 3.5.0)
##  graphics                    * 3.5.0    2017-05-24 local         
##  grDevices                   * 3.5.0    2017-05-24 local         
##  grid                          3.5.0    2017-05-24 local         
##  here                        * 0.1      2017-05-28 CRAN (R 3.5.0)
##  hms                           0.3      2016-11-22 CRAN (R 3.5.0)
##  htmltools                     0.3.6    2017-04-28 CRAN (R 3.5.0)
##  httpuv                        1.3.5    2017-07-04 CRAN (R 3.5.0)
##  httr                          1.3.1    2017-08-20 CRAN (R 3.5.0)
##  interactiveDisplayBase        1.15.0   2017-08-12 Bioconductor  
##  IRanges                     * 2.11.19  2017-10-18 Bioconductor  
##  jsonlite                      1.5      2017-06-01 CRAN (R 3.5.0)
##  knitr                         1.17     2017-08-10 CRAN (R 3.5.0)
##  lattice                       0.20-35  2017-03-25 CRAN (R 3.5.0)
##  lazyeval                      0.2.0    2016-06-12 CRAN (R 3.5.0)
##  magrittr                    * 1.5      2014-11-22 CRAN (R 3.5.0)
##  Matrix                        1.2-11   2017-08-16 CRAN (R 3.5.0)
##  matrixStats                 * 0.52.2   2017-04-14 CRAN (R 3.5.0)
##  memoise                       1.1.0    2017-04-21 CRAN (R 3.5.0)
##  methods                     * 3.5.0    2017-05-24 local         
##  mime                          0.5      2016-07-07 CRAN (R 3.5.0)
##  parallel                    * 3.5.0    2017-05-24 local         
##  pkgconfig                     2.0.1    2017-03-21 CRAN (R 3.5.0)
##  prettyunits                   1.0.2    2015-07-13 CRAN (R 3.5.0)
##  progress                      1.1.2    2016-12-14 CRAN (R 3.5.0)
##  ProtGenerics                  1.9.1    2017-10-08 Bioconductor  
##  R6                            2.2.2    2017-06-17 CRAN (R 3.5.0)
##  rappdirs                      0.3.1    2016-03-28 CRAN (R 3.5.0)
##  Rcpp                          0.12.13  2017-09-28 CRAN (R 3.5.0)
##  RCurl                         1.95-4.8 2016-03-01 CRAN (R 3.5.0)
##  readr                       * 1.1.1    2017-05-16 CRAN (R 3.5.0)
##  rjson                         0.2.15   2014-11-03 CRAN (R 3.5.0)
##  rlang                         0.1.2    2017-08-09 CRAN (R 3.5.0)
##  rmarkdown                   * 1.6      2017-06-15 CRAN (R 3.5.0)
##  roxygen2                      6.0.1    2017-02-06 CRAN (R 3.5.0)
##  rprojroot                     1.2      2017-01-16 CRAN (R 3.5.0)
##  Rsamtools                     1.29.1   2017-10-08 Bioconductor  
##  RSQLite                       2.0      2017-06-19 CRAN (R 3.5.0)
##  rstudioapi                    0.7      2017-09-07 CRAN (R 3.5.0)
##  rtracklayer                 * 1.37.3   2017-10-08 Bioconductor  
##  S4Vectors                   * 0.15.14  2017-10-18 Bioconductor  
##  shiny                         1.0.5    2017-08-23 CRAN (R 3.5.0)
##  stats                       * 3.5.0    2017-05-24 local         
##  stats4                      * 3.5.0    2017-05-24 local         
##  stringi                       1.1.5    2017-04-07 CRAN (R 3.5.0)
##  stringr                       1.2.0    2017-02-18 CRAN (R 3.5.0)
##  SummarizedExperiment        * 1.7.10   2017-10-08 Bioconductor  
##  testthat                    * 1.0.2    2016-04-23 CRAN (R 3.5.0)
##  tibble                        1.3.4    2017-08-22 CRAN (R 3.5.0)
##  tools                         3.5.0    2017-05-24 local         
##  tximeta                     * 0.0.5    <NA>       Bioconductor  
##  tximport                      1.5.1    2017-10-08 Bioconductor  
##  utils                       * 3.5.0    2017-05-24 local         
##  withr                         2.0.0    2017-07-28 CRAN (R 3.5.0)
##  XML                           3.98-1.9 2017-06-19 CRAN (R 3.5.0)
##  xml2                          1.1.1    2017-01-24 CRAN (R 3.5.0)
##  xtable                        1.8-2    2016-02-05 CRAN (R 3.5.0)
##  XVector                     * 0.17.1   2017-10-08 Bioconductor  
##  yaml                          2.1.14   2016-11-12 CRAN (R 3.5.0)
##  zlibbioc                      1.23.0   2017-05-24 Bioconductor
```
