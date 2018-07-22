---
title: "tximeta: transcript quantification import with automatic metadata"
author: "Michael Love, Rob Patro, Charlotte Soneson, Peter Hickey"
date: "07/22/2018"
output: 
  rmarkdown::html_document:
    highlight: tango
abstract: >
  `tximeta` performs numerous annotation and metadata gathering tasks on
  behalf of users during the import of transcript quantifications from
  *Salmon* or *Sailfish* into R/Bioconductor. Metadata and transcript
  ranges are added automatically, facilitating combining multiple
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
##                                                                                  files
## 1 /Users/love/bin/R/library/tximportData/extdata/salmon_dm/SRR1197474_cdna/quant.sf.gz
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

However, to avoid downloading remote GTF files during this vignette,
we will point to a GTF file saved locally (in the *tximportData*
package). We link the transcriptome of the *Salmon* index to its
locally saved GTF. The standard recommended usage of `tximeta` would
be the code chunk above, or to specify a remote GTF source, not a
local one. This following code is therefore not recommended for a
typically workflow, but is particular to the vignette code.


```r
dir <- system.file("extdata", package="tximeta")
indexDir <- file.path(dir, "Drosophila_melanogaster.BDGP6.cdna.v92_salmon_0.10.2")
fastaFTP <- "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
dir2 <- system.file("extdata/salmon_dm", package="tximportData")
gtfPath <- file.path(dir2,"Drosophila_melanogaster.BDGP6.92.gtf.gz")
suppressPackageStartupMessages(library(tximeta))
makeLinkedTxome(indexDir=indexDir,
                source="Ensembl",
                organism="Drosophila melanogaster",
                release="92",
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
## [ Ensembl - Drosophila melanogaster - release 92 ]
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
index before, it will use a cached version of the metadata and
transcript ranges. 

Note the warning above that 9 of the transcripts are missing from the
GTF file and so are dropped from the final output. This is a problem
coming from the annotation source, and not easily avoided by
`tximeta`. 

We plan to create and maintain a large table of signatures for as many
sources, organisms, releases of transcriptomes as possible. We are
also developing support for *linked transcriptomes*, where one or
more sources for transcript sequences have been combined or
filtered. See the `linkedTxome` vignette in this package for a
demonstration. (The *makeLinkedTxome* function was used above to avoid
downloading the GTF during the vignette building process.)

# Examining SummarizedExperiment output

We, of course, have our coldata from before. Note that we've removed `files`.


```r
suppressPackageStartupMessages(library(SummarizedExperiment))
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
##                 gene_id   gene_name  entrezid   gene_biotype
##             <character> <character> <integer>    <character>
## FBgn0000008 FBgn0000008           a        NA protein_coding
## FBgn0000014 FBgn0000014       abd-A        NA protein_coding
## FBgn0000015 FBgn0000015       Abd-B        NA protein_coding
## FBgn0000017 FBgn0000017         Abl        NA protein_coding
## FBgn0000018 FBgn0000018         abo        NA protein_coding
## ...                 ...         ...       ...            ...
## FBgn0285958 FBgn0285958        Fuca        NA protein_coding
## FBgn0285962 FBgn0285962     CG46338        NA protein_coding
## FBgn0285963 FBgn0285963     CG46339        NA protein_coding
## FBgn0285970 FBgn0285970     CG32500        NA protein_coding
## FBgn0285971 FBgn0285971         prg        NA protein_coding
##             seq_coord_system      symbol       REFSEQ
##                    <integer> <character>  <character>
## FBgn0000008               NA           a NM_001014543
## FBgn0000014               NA       abd-A NM_001170161
## FBgn0000015               NA       Abd-B NM_001275719
## FBgn0000017               NA         Abl NM_001104153
## FBgn0000018               NA         abo    NM_080045
## ...                      ...         ...          ...
## FBgn0285958               NA        Fuca NM_001316434
## FBgn0285962               NA     CG46338 NM_001273908
## FBgn0285963               NA     CG46339 NM_001104469
## FBgn0285970               NA     CG32500    NM_167766
## FBgn0285971               NA         prg NM_001273324
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
##  $ release       : chr "92"
##  $ genome        : chr "BDGP6"
##  $ fasta         :List of 1
##   ..$ : chr "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
##  $ gtf           : chr "/Users/love/bin/R/library/tximportData/extdata/salmon_dm/Drosophila_melanogaster.BDGP6.92.gtf.gz"
```

```r
str(metadata(se)$tximetaInfo)
```

```
## List of 2
##  $ version   :Classes 'package_version', 'numeric_version'  hidden list of 1
##   ..$ : int [1:3] 0 99 6
##  $ importTime: POSIXct[1:1], format: "2018-07-22 07:42:43"
```

```r
str(metadata(se)$txdbInfo)
```

```
##  Named chr [1:12] "EnsDb" "Ensembl Gene ID" "ensembldb" ...
##  - attr(*, "names")= chr [1:12] "Db type" "Type of Gene ID" "Supporting package" "Db created by" ...
```

# Quantification files with an unknown transcriptome

`tximeta` automatically imports relevant metadata when the
transcriptome matches a known source, but also facilitates the
linking of transcriptomes used as for a *Salmon* index with relevant
public sources. The linking is important in the case that the
transcript sequence no longer matches a known source (combined or
filtered FASTA files), or if the source is not known to
`tximeta`. Below we demonstrate how to make a *linkedTxome* and how
to share and load a *linkedTxome*.

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

# Linked transcriptome for reproducible analysis

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
non-coding FASTA files *Drosophila melanogaster*, release 92. The
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

Typically you would not use `system.file` to find this directory, but
simply define `indexDir` to be the path of the *Salmon* directory on
your machine. Here we use `system.file` because we have included parts
of a *Salmon* index directory in the *tximeta* package itself for
demonstration of functionality in this vignette.


```r
dir <- system.file("extdata", package="tximeta")
indexDir <- file.path(dir, "Drosophila_melanogaster.BDGP6.v92_salmon_0.10.2")
```

Now we provide the location of the FASTA files and the GTF file for
this transcriptome. The recommended usage of `tximeta` would be to
specify a remote GTF source, as seen in the commented-out line below:


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
                release="92", genome="BDGP6",
                fasta=fastaFTP, gtf=gtfPath,
                jsonFile=jsonFile)
```

```
## writing linkedTxome to /var/folders/kd/xns8y9c51sz668rrjn0bvlqh0000gn/T//Rtmp7XCj6Q/Drosophila_melanogaster.BDGP6.v92_salmon_0.10.2.json
```

```
## saving linkedTxome in bfc
```

After running `makeLinkedTxome`, the connection between this *Salmon*
index (and its signature) with the sources is saved for persistent
usage.

With use of `tximeta` and a *linkedTxome* -- as with `tximeta` on a
known, un-filtered, un-combined transcriptome -- the software
figures out if the remote GTF has been accessed and compiled into a
*TxDb* before, and on future calls, it will simply load the
pre-computed metadata and transcript ranges.

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
## [ Ensembl - Drosophila melanogaster - release 92 ]
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

The following code removes the entire table with information about the
*linkedTxomes*. This is just for demonstration, so that we can show
how to load a JSON file below.

**Note:** Running this code will clear any information about
*linkedTxomes*. Don't run this unless you really want to clear this
table!


```r
library(BiocFileCache)
if (interactive()) {
  bfcloc <- getTximetaBFC()
} else {
  bfcloc <- tempdir()
}
bfc <- BiocFileCache(bfcloc)
bfcinfo(bfc)
```

```
## # A tibble: 3 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## 2 BFC41 genco… 2018-07-19… 2018-07-19… /User… rela… bf40…               NA
## 3 BFC46 linke… 2018-07-22… 2018-07-22… /User… rela… 9c3d…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

```r
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 2 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## 2 BFC41 genco… 2018-07-19… 2018-07-19… /User… rela… bf40…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

# Loading *linkedTxome* JSON files

If a collaborator or the Suppmentary Files for a publication shares a
`linkedTxome` JSON file, we can likewise use `tximeta` to
automatically assemble the relevant metadata and transcript
ranges. This implies that the other person has used `tximeta` with the
function `makeLinkedTxome` demonstrated above, pointing to their
*Salmon* index and to the FASTA and GTF source(s).

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
## [ Ensembl - Drosophila melanogaster - release 92 ]
## loading existing EnsDb created: 2018-07-14 14:45:15
## generating transcript ranges
```

```
## Warning in checkTxi2Txps(txi, txps): missing some transcripts!
##  9 out of 33681 are missing from the GTF and dropped from SummarizedExperiment output
```

# Clear *linkedTxomes*

Finally, we clear the *linkedTxomes* table again so that the above
examples will work. This is just for the vignette code and not part of
a typical workflow.

**Note:** Running this code will clear any information about
*linkedTxomes*. Don't run this unless you really want to clear this
table!


```r
if (interactive()) {
  bfcloc <- getTximetaBFC()
} else {
  bfcloc <- tempdir()
}
bfc <- BiocFileCache(bfcloc)
bfcinfo(bfc)
```

```
## # A tibble: 3 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## 2 BFC41 genco… 2018-07-19… 2018-07-19… /User… rela… bf40…               NA
## 3 BFC47 linke… 2018-07-22… 2018-07-22… /User… rela… 9c33…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
```

```r
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcinfo(bfc)
```

```
## # A tibble: 2 x 10
##   rid   rname  create_time access_time rpath  rtype fpath last_modified_t…
##   <chr> <chr>  <chr>       <chr>       <chr>  <chr> <chr>            <dbl>
## 1 BFC28 Droso… 2018-07-14… 2018-07-14… /User… rela… 31e7…               NA
## 2 BFC41 genco… 2018-07-19… 2018-07-19… /User… rela… bf40…               NA
## # ... with 2 more variables: etag <chr>, expires <dbl>
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
  `release` (also here we ignored something like "type", e.g. CHR
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
##  version  R Under development (unstable) (2018-05-14 r74725)
##  system   x86_64, darwin15.6.0                              
##  ui       X11                                               
##  language (EN)                                              
##  collate  en_US.UTF-8                                       
##  tz       America/New_York                                  
##  date     2018-07-22
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version   date       source        
##  acepack                1.4.1     2016-10-29 CRAN (R 3.6.0)
##  annotate               1.59.0    2018-05-15 Bioconductor  
##  AnnotationDbi        * 1.43.1    2018-05-15 Bioconductor  
##  AnnotationFilter       1.5.2     2018-06-13 Bioconductor  
##  assertthat             0.2.0     2017-04-11 CRAN (R 3.6.0)
##  backports              1.1.2     2017-12-13 CRAN (R 3.6.0)
##  base                 * 3.6.0     2018-05-15 local         
##  base64enc              0.1-3     2015-07-28 CRAN (R 3.6.0)
##  bindr                  0.1.1     2018-03-13 CRAN (R 3.6.0)
##  bindrcpp             * 0.2.2     2018-03-29 CRAN (R 3.6.0)
##  Biobase              * 2.41.1    2018-06-17 Bioconductor  
##  BiocFileCache        * 1.5.4     2018-07-06 Bioconductor  
##  BiocGenerics         * 0.27.1    2018-06-17 Bioconductor  
##  BiocInstaller        * 1.31.1    2018-05-15 Bioconductor  
##  BiocParallel         * 1.15.6    2018-06-28 Bioconductor  
##  biomaRt                2.37.3    2018-06-29 Bioconductor  
##  Biostrings             2.49.0    2018-05-22 Bioconductor  
##  bit                    1.1-14    2018-05-29 CRAN (R 3.6.0)
##  bit64                  0.9-7     2017-05-08 CRAN (R 3.6.0)
##  bitops                 1.0-6     2013-08-17 CRAN (R 3.6.0)
##  blob                   1.1.1     2018-03-25 CRAN (R 3.6.0)
##  checkmate              1.8.5     2017-10-24 CRAN (R 3.6.0)
##  cli                    1.0.0     2017-11-05 CRAN (R 3.6.0)
##  cluster                2.0.7-1   2018-04-13 CRAN (R 3.6.0)
##  colorspace             1.3-2     2016-12-14 CRAN (R 3.6.0)
##  compiler               3.6.0     2018-05-15 local         
##  crayon                 1.3.4     2017-09-16 CRAN (R 3.6.0)
##  curl                   3.2       2018-03-28 CRAN (R 3.6.0)
##  data.table             1.11.4    2018-05-27 CRAN (R 3.6.0)
##  datasets             * 3.6.0     2018-05-15 local         
##  DBI                    1.0.0     2018-05-02 CRAN (R 3.6.0)
##  dbplyr               * 1.2.1     2018-02-19 CRAN (R 3.6.0)
##  DelayedArray         * 0.7.19    2018-07-03 Bioconductor  
##  DESeq2               * 1.21.7    2018-07-03 Bioconductor  
##  devtools             * 1.13.6    2018-06-27 CRAN (R 3.6.0)
##  digest                 0.6.15    2018-01-28 CRAN (R 3.6.0)
##  dplyr                  0.7.6     2018-06-29 CRAN (R 3.6.0)
##  ensembldb              2.5.3     2018-07-05 Bioconductor  
##  evaluate               0.10.1    2017-06-24 CRAN (R 3.6.0)
##  foreign                0.8-70    2017-11-28 CRAN (R 3.6.0)
##  Formula                1.2-3     2018-05-03 CRAN (R 3.6.0)
##  genefilter             1.63.0    2018-05-15 Bioconductor  
##  geneplotter            1.59.0    2018-05-15 Bioconductor  
##  GenomeInfoDb         * 1.17.1    2018-05-15 Bioconductor  
##  GenomeInfoDbData       1.1.0     2018-05-15 Bioconductor  
##  GenomicAlignments      1.17.2    2018-06-13 Bioconductor  
##  GenomicFeatures        1.33.0    2018-06-13 Bioconductor  
##  GenomicRanges        * 1.33.6    2018-06-13 Bioconductor  
##  ggplot2                3.0.0     2018-07-03 CRAN (R 3.6.0)
##  glue                   1.2.0     2017-10-29 CRAN (R 3.6.0)
##  graphics             * 3.6.0     2018-05-15 local         
##  grDevices            * 3.6.0     2018-05-15 local         
##  grid                   3.6.0     2018-05-15 local         
##  gridExtra              2.3       2017-09-09 CRAN (R 3.6.0)
##  gtable                 0.2.0     2016-02-26 CRAN (R 3.6.0)
##  Hmisc                  4.1-1     2018-01-03 CRAN (R 3.6.0)
##  hms                    0.4.2     2018-03-10 CRAN (R 3.6.0)
##  htmlTable              1.12      2018-05-26 CRAN (R 3.6.0)
##  htmltools              0.3.6     2017-04-28 CRAN (R 3.6.0)
##  htmlwidgets            1.2       2018-04-19 CRAN (R 3.6.0)
##  httr                   1.3.1     2017-08-20 CRAN (R 3.6.0)
##  IRanges              * 2.15.14   2018-06-13 Bioconductor  
##  jsonlite               1.5       2017-06-01 CRAN (R 3.6.0)
##  knitr                  1.20      2018-02-20 CRAN (R 3.6.0)
##  lattice                0.20-35   2017-03-25 CRAN (R 3.6.0)
##  latticeExtra           0.6-28    2016-02-09 CRAN (R 3.6.0)
##  lazyeval               0.2.1     2017-10-29 CRAN (R 3.6.0)
##  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.6.0)
##  magrittr               1.5       2014-11-22 CRAN (R 3.6.0)
##  Matrix                 1.2-14    2018-04-13 CRAN (R 3.6.0)
##  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.6.0)
##  memoise                1.1.0     2017-04-21 CRAN (R 3.6.0)
##  methods              * 3.6.0     2018-05-15 local         
##  munsell                0.5.0     2018-06-12 CRAN (R 3.6.0)
##  nnet                   7.3-12    2016-02-02 CRAN (R 3.6.0)
##  org.Dm.eg.db         * 3.6.0     2018-07-14 Bioconductor  
##  parallel             * 3.6.0     2018-05-15 local         
##  pillar                 1.2.3     2018-05-25 CRAN (R 3.6.0)
##  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.6.0)
##  plyr                   1.8.4     2016-06-08 CRAN (R 3.6.0)
##  prettyunits            1.0.2     2015-07-13 CRAN (R 3.6.0)
##  progress               1.2.0     2018-06-14 CRAN (R 3.6.0)
##  ProtGenerics           1.13.0    2018-06-13 Bioconductor  
##  purrr                  0.2.5     2018-05-29 CRAN (R 3.6.0)
##  R6                     2.2.2     2017-06-17 CRAN (R 3.6.0)
##  rappdirs               0.3.1     2016-03-28 CRAN (R 3.6.0)
##  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.6.0)
##  Rcpp                   0.12.16   2018-03-13 CRAN (R 3.6.0)
##  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.6.0)
##  readr                  1.1.1     2017-05-16 CRAN (R 3.6.0)
##  rlang                  0.2.0     2018-02-20 CRAN (R 3.6.0)
##  rmarkdown            * 1.10      2018-06-11 CRAN (R 3.6.0)
##  rpart                  4.1-13    2018-02-23 CRAN (R 3.6.0)
##  rprojroot              1.3-2     2018-01-03 CRAN (R 3.6.0)
##  Rsamtools              1.33.2    2018-07-03 Bioconductor  
##  RSQLite                2.1.1     2018-05-06 CRAN (R 3.6.0)
##  rstudioapi             0.7       2017-09-07 CRAN (R 3.6.0)
##  rtracklayer            1.41.3    2018-06-13 Bioconductor  
##  S4Vectors            * 0.19.17   2018-06-28 Bioconductor  
##  scales                 0.5.0     2017-08-24 CRAN (R 3.6.0)
##  splines                3.6.0     2018-05-15 local         
##  stats                * 3.6.0     2018-05-15 local         
##  stats4               * 3.6.0     2018-05-15 local         
##  stringi                1.2.2     2018-05-02 CRAN (R 3.6.0)
##  stringr                1.3.1     2018-05-10 CRAN (R 3.6.0)
##  SummarizedExperiment * 1.11.5    2018-06-13 Bioconductor  
##  survival               2.42-4    2018-06-30 CRAN (R 3.6.0)
##  testthat             * 2.0.0     2017-12-13 CRAN (R 3.6.0)
##  tibble                 1.4.2     2018-01-22 CRAN (R 3.6.0)
##  tidyselect             0.2.4     2018-02-26 CRAN (R 3.6.0)
##  tools                  3.6.0     2018-05-15 local         
##  tximeta              * 0.99.6    2018-07-22 Bioconductor  
##  tximport               1.9.8     2018-07-11 cran (@1.9.8) 
##  utf8                   1.1.4     2018-05-24 CRAN (R 3.6.0)
##  utils                * 3.6.0     2018-05-15 local         
##  withr                  2.1.2     2018-03-15 CRAN (R 3.6.0)
##  XML                    3.98-1.11 2018-04-16 CRAN (R 3.6.0)
##  xtable                 1.8-2     2016-02-05 CRAN (R 3.6.0)
##  XVector                0.21.3    2018-06-23 Bioconductor  
##  yaml                   2.1.19    2018-05-01 CRAN (R 3.6.0)
##  zlibbioc               1.27.0    2018-05-15 Bioconductor
```
