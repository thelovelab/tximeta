# Make and load linked transcriptomes (linked GTF and FASTA)

`makeLinkedTxome()` reads the digest associated with a salmon index at
`indexDir`, and persistently links it to metadata (alternatively the
`digest` string itself and an `indexName` can be provided). Linked
metadata includes key information about the transcriptome, including the
`source`, `organism`, `release`, and `genome` (these are custom
character strings), as well as the locations (e.g. local, HTTP, or FTP)
for one or more `fasta` files and one `gtf` file. `loadLinkedTxome()`
loads this information from a JSON file. See *Details*.

## Usage

``` r
makeLinkedTxome(
  digest = NULL,
  indexName,
  indexDir = NULL,
  source,
  organism,
  release,
  genome,
  fasta,
  gtf,
  write = TRUE,
  jsonFile
)

loadLinkedTxome(jsonFile)
```

## Arguments

- digest:

  the full digest as character string, (this or `indexDir` is required,
  only one should be specified)

- indexName:

  a name for the `index` when storing the linkedTxome, required if
  providing the `digest` string, suggest using the basename of the FASTA
  file and the software used, e.g. "gencode.vXX_salmon-0.XX.Y"

- indexDir:

  the local path to the salmon index (this or `digest` is required, only
  one should be specified)

- source:

  the source of transcriptome (e.g. "de-novo"). Note: if you specify
  "GENCODE" or "Ensembl", this will trigger behavior by tximeta that may
  not be desired: e.g. attempts to download canonical transcriptome data
  from AnnotationHub (unless useHub=FALSE when running tximeta) and
  parsing of Ensembl GTF using ensembldb (which may fail if the GTF file
  has been modified). For transcriptomes that are defined by local GTF
  files, it is recommended to use the terms "LocalGENCODE" or
  "LocalEnsembl". Setting "LocalEnsembl" will also strip version numbers
  from the FASTA transcript IDs to enable matching with the Ensembl GTF.

- organism:

  organism (e.g. "Homo sapiens")

- release:

  release number (e.g. "27")

- genome:

  genome (e.g. "GRCh38", or "none")

- fasta:

  location(s) for the FASTA transcript sequences (of which the
  transcripts used to build the index is equal or a subset). This can be
  a local path, or an HTTP or FTP URL

- gtf:

  location for the GTF/GFF file (of which the transcripts used to build
  the index is equal or a subset). This can be a local path, or an HTTP
  or FTP URL While the `fasta` argument can take a vector of length
  greater than one (more than one FASTA file containing transcripts used
  in indexing), the `gtf` argument has to be a single GTF/GFF file. This
  can also be a serialized GRanges object (location of a .rds file)
  imported with rtracklayer. If transcripts were added to a standard set
  of reference transcripts (e.g. fusion genes, or pathogen transcripts),
  it is recommended that the tximeta user would manually add these to
  the GTF/GFF file, and post the modified GTF/GFF publicly, such as on
  Zenodo. This enables consistent annotation and downstream annotation
  tasks, such as by
  [`summarizeToGene()`](https://thelovelab.github.io/tximeta/reference/summarizeToGene.md).

- write:

  logical, should a JSON file be written out which documents the
  transcriptome digest and metadata? (default is TRUE)

- jsonFile:

  the path to the json file for the linkedTxome

## Value

nothing, the function is run for its side effects

## Details

`makeLinkedTxome()` links the information about the transcriptome used
for quantification in two ways:

1.  the function will store a record in tximeta's cache such that future
    import of quantification data will automatically access and parse
    the GTF as if the transcriptome were one of those automatically
    detected by tximeta. Then all features of tximeta (e.g.
    summarization to gene, programmatic adding of IDs or metadata) will
    be available;

2.  it will by default write out a JSON file that can be shared, or
    posted online, and which can be read by `loadLinkedTxome()` which
    will store the information in tximeta's cache. This should make the
    full quantification-import pipeline computationally reproducible /
    auditable even for transcriptomes which differ from those provided
    by references (GENCODE, Ensembl, RefSeq).

For further details please see the "Linked transcriptomes" section of
the tximeta vignette.

This function can be used in combination with
[`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
and oarfish data from
[`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md),
when multiple reference transcript sets have been indexed. See also
[`makeLinkedTxpData()`](https://thelovelab.github.io/tximeta/reference/linkedTxpData.md).

## Examples

``` r

# point to a salmon quantification file with an additional artificial transcript
dir <- system.file("extdata/salmon_dm", package="tximportData")
file <- file.path(dir, "SRR1197474.plus", "quant.sf")
coldata <- data.frame(files=file, names="SRR1197474", sample="1",
                      stringsAsFactors=FALSE)

# now point to the salmon index itself to create a linkedTxome
# as the index will not match a known txome
indexDir <- file.path(dir, "Dm.BDGP6.22.98.plus_salmon-0.14.1")

# point to the source FASTA and GTF:
baseFTP <- "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/"
fastaFTP <- c(
  paste0(baseFTP,
    c("cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
      "ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")),
  "extra_transcript.fa.gz"
)
gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.plus.gtf.gz")

# now create a linkedTxome, linking the salmon index to its FASTA and GTF sources
makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster",
                release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)
#> reading digest from indexDir: .../Dm.BDGP6.22.98.plus_salmon-0.14.1
#> saving linkedTxome in bfc

# to clear the entire linkedTxome table
# (don't run unless you want to clear this table!)
# bfcloc <- getTximetaBFC()
# bfc <- BiocFileCache(bfcloc)
# bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
```
