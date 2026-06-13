# Import transcript quantification with metadata

`tximeta` leverages the digest of the reference transcripts that were
indexed in order to identify metadata from the output of quantification
tools. A computed digest (a hash value) can be used to uniquely identify
the collection of reference sequences, and associate the dataset with
other useful metadata. After identification, tximeta uses a number of
core Bioconductor packages (GenomicFeatures, ensembldb, AnnotationHub,
Seqinfo, BiocFileCache) to automatically populate metadata for the user.

## Usage

``` r
tximeta(
  coldata,
  type = NULL,
  txOut = TRUE,
  skipMeta = FALSE,
  skipSeqinfo = FALSE,
  useHub = TRUE,
  markDuplicateTxps = FALSE,
  cleanDuplicateTxps = FALSE,
  customMetaInfo = NULL,
  skipFtp = FALSE,
  gencode_gtf_prefix = NULL,
  ...
)
```

## Arguments

- coldata:

  a data.frame with at least two columns (others will propogate to
  object):

  - `files` - character, paths of quantification files

  - `names` - character, sample names if `coldata` is a vector, it is
    assumed to be the paths of quantification files and unique sample
    names are created

- type:

  what quantifier was used, see
  [`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html)

- txOut:

  whether to output transcript-level data. `tximeta` is designed to have
  transcript-level output with salmon, so default is `TRUE`, and it's
  recommended to use `summarizeToGene` following `tximeta` for
  gene-level summarization. For an alevin file, `tximeta` will import
  the gene level counts ignoring this argument (alevin produces only
  gene-level quantification).

- skipMeta:

  whether to skip metadata generation (e.g. to avoid errors if not
  connected to internet). This calls `tximport` directly and so either
  `txOut=TRUE` or `tx2gene` should be specified.

- skipSeqinfo:

  whether to skip the addition of Seqinfo, which requires an internet
  connection to download the relevant chromosome information table from
  UCSC

- useHub:

  whether to first attempt to download a TxDb/EnsDb object from
  AnnotationHub, rather than creating from a GTF file from FTP (default
  is TRUE). If FALSE, it will force `tximeta` to download and parse the
  GTF

- markDuplicateTxps:

  whether to mark the status (`hasDuplicate`) and names of duplicate
  transcripts (`duplicates`) in the rowData of the SummarizedExperiment
  output. Subsequent summarization to gene level will keep track of the
  number of transcripts sets per gene (`numDupSets`). see Details

- cleanDuplicateTxps:

  whether to "clean" duplicate transcripts (exact sequence duplicates)
  by replacing, when possible, transcript names that do not appear in
  the GTF with those that do appear in the GTF. see Details

- customMetaInfo:

  the relative path to a custom metadata information JSON file, relative
  to the paths in `files` of `coldata`. For example,
  `customMetaInfo="meta_info.json"` would indicate that in the same
  directory as the quantification files in `files`, there are custom
  metadata information JSON files. These should contain the SHA-256 hash
  of the reference transcripts with the `index_seq_hash` tag (see
  details in vignette).

- skipFtp:

  whether to avoid `ftp://` in case of firewall, default is FALSE

- gencode_gtf_prefix:

  for GENCODE transcriptomes, optionally specify a non-default GTF file
  by providing the prefix that appears between `gencode.vXX.` and
  `.annotation.gtf.gz` in the GENCODE FTP filename. For example,
  `"primary_assembly"`, `"basic"`, `"chr_patch_hapl_scaff"`,
  `"primary_assembly.basic"`, `"chr_patch_hapl_scaff.basic"`. Has no
  effect for non-GENCODE transcriptomes.

- ...:

  arguments passed to `tximport`

## Value

a SummarizedExperiment with metadata on the `rowRanges`. (if the hashed
digest in the salmon or Sailfish index does not match any known
transcriptomes, or any locally saved `linkedTxome`, `tximeta` will just
return a non-ranged SummarizedExperiment)

## Details

Most of the code in tximeta works to add metadata and transcript ranges
when the quantification was performed with salmon or related tools.
However, tximeta can be used with any quantification type that is
supported by
[`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html),
where it will return an non-ranged SummarizedExperiment. For other
quantification tools see also the `customMetaInfo` argument below. This
behavior can also be triggered with `skipMeta=TRUE`.

tximeta performs a lookup of the digest (or hash value) of the index
stored in an auxilary information directory of the quantification tool's
output against a database of known transcriptomes, which is stored
within the tximeta package (`extdata/hashtable.csv`) and is continually
updated to match Ensembl and GENCODE releases, with updates pushed to
Bioconductor current release branch. In addition, tximeta performs a
lookup of the digest against a locally stored table of linkedTxome
references, see
[`makeLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md).
If tximeta detects a match in either source, it will automatically
populate the transcript locations, the transcriptome release, the genome
with correct chromosome lengths, and connect the SE object to locally
cached derived metadata. tximeta also facilitates automatic
summarization of transcript-level quantifications to the gene-level via
``` summarizeToGene`` without the need to manually build the correct  ```tx2gene\`
table for the reference used for indexing.

tximeta on the first run will ask where the
[`BiocFileCache::BiocFileCache()`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
location for this package (*tximeta*) should be kept, either using a
default location or a temporary directory. At any point, the user can
specify a location using
[`setTximetaBFC()`](https://thelovelab.github.io/tximeta/reference/getTximetaBFC.md)
and this choice will be saved for future sessions. Multiple users can
point to the same BiocFileCache, such that transcript databases (TxDb or
EnsDb) associated with certain salmon indices and linkedTxomes can be
accessed by different users without additional effort or time spent
downloading and building the relevant TxDb / EnsDb. In order to allow
that multiple users can read and write to the same location, one should
set the BiocFileCache directory to have group write permissions (g+w).
Note that, if the TxDb or EnsDb is present in AnnotationHub, tximeta
will use this object instead of downloading and building a TxDb/EnsDb
from GTF (to disable this set `useHub=FALSE`).

Regarding `markDuplicateTxps` and `cleanDuplicateTxps`: these functions
may be used when salmon has been run in default mode (where it reduces
exact sequence duplicates in the index by replacing them with a single
representative). These won't make sense (and the latter will give an
error) when `--keepDuplicates` has been used by salmon during indexing.
They help by either adding metadata ("marking") regarding transcript
duplicates to the rowData, or changing names of transcripts ("cleaning")
where possible to provide better alignment with GTF.

## Examples

``` r

# point to a salmon quantification file:
dir <- system.file("extdata/salmon_dm", package="tximportData")
files <- file.path(dir, "SRR1197474", "quant.sf") 
coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)

# normally we would just run the following which would download the appropriate metadata
# se <- tximeta(coldata)

# for this example, we instead point to a local path where the GTF can be found
# by making a linkedTxome:
indexDir <- file.path(dir, "Dm.BDGP6.22.98_salmon-0.14.1")
dmFTP <- "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/"
fastaFTP <- paste0(
  dmFTP,
  c("cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
    "ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")
)
gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster",
                release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#> reading digest from indexDir: .../Dm.BDGP6.22.98_salmon-0.14.1
#> NOTE: this digest matches one in the pre-computed digest table
#> linkedTxome metadata was same as already in bfc
se <- tximeta(coldata)
#> importing salmon quantification files
#> reading in files with read.delim (install 'readr' package for speed up)
#> 1 
#> 
#> found matching linkedTxome:
#> [ LocalEnsembl - Drosophila melanogaster - release 98 ]
#> loading existing TxDb created: 2026-06-13 17:45:19
#> loading existing transcript ranges created: 2026-06-13 17:45:19
#> Warning: 
#> 
#> Warning: the annotation is missing some transcripts that were quantified.
#> 5 out of 33706 txps were missing from GTF/GFF but were in the indexed FASTA
#> (e.g. this can occur with transcripts located on haplotype chromosomes).
#> In order to build a ranged SummarizedExperiment, these txps were removed.
#> To keep these txps, and to skip adding ranges, use skipMeta=TRUE
#> 
#> Example missing txps: [FBtr0307759, FBtr0084079, FBtr0084080, ...]

# to clear the entire linkedTxome table
# (don't run unless you want to clear this table!)
# bfcloc <- getTximetaBFC()
# bfc <- BiocFileCache(bfcloc)
# bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
```
