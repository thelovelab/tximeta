---
title: "tximeta: Working with derived transcriptomes"
author: "Michael Love, Rob Patro"
date: "12/02/2017"
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
download.file("https://github.com/mikelove/crohns/archive/master.zip", dest)
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
`index` step. By default, `derivedTxome` will write out a JSON file
which can be shared with others, linking the signature of the index
with the other metadata, including FASTA and GTF sources. It will
write out to a file with the same name as the `indexDir`, but with
a `.json` extension added.


```r
makeDerivedTxome(
  indexDir="Homo_sapiens.GRCh38.cdna.std.chroms_salmon_0.8.2",
  source="Ensembl",
  organism="Homo sapiens",
  version="90",
  genome="GRCh38",
  fasta="ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
  gtf="ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz"
)
```

```
## writing derivedTxome to Homo_sapiens.GRCh38.cdna.std.chroms_salmon_0.8.2.json
```

```
## saving derivedTxome in bfc (first time)
```

Now it works!


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
## loading existing TxDb created: 2017-12-03 01:36:31
## generating transcript ranges
## fetching genome info
```

Remove to make it not work from the beginning


```r
library(BiocFileCache)
bfc <- BiocFileCache(".")
bfcinfo(bfc)
```

```
## # A tibble: 2 x 8
##     rid                         rname         create_time
##   <chr>                         <chr>               <chr>
## 1 BFC11 Homo_sapiens.GRCh38.90.gtf.gz 2017-12-03 01:36:31
## 2 BFC13                derivedTxomeDF 2017-12-03 01:43:26
## # ... with 5 more variables: access_time <chr>, rpath <chr>, rtype <chr>,
## #   fpath <chr>, last_modified_time <chr>
```

```r
bfcremove(bfc, bfcquery(bfc, "derivedTxomeDF")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 1 x 8
##     rid                         rname         create_time
##   <chr>                         <chr>               <chr>
## 1 BFC11 Homo_sapiens.GRCh38.90.gtf.gz 2017-12-03 01:36:31
## # ... with 5 more variables: access_time <chr>, rpath <chr>, rtype <chr>,
## #   fpath <chr>, last_modified_time <chr>
```
