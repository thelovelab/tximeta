# Make linked transcript data (linked GRanges)

`linkedTxpData` allows the user to save relevant *GRanges* transcript
data for identifying and updating transcript metadata in a persistent
manner across R sessions. It can be used in combination with
[`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
and
[`updateMetadata()`](https://thelovelab.github.io/tximeta/reference/updateMetadata.md).
This is a lightweight version of `linkedTxome` (see
[`makeLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md)),
which requires specifying a GTF file for building a *TxDb* and
optionally a FASTA file for sequence retrieval.)

## Usage

``` r
makeLinkedTxpData(
  digest,
  digestType = "sha256",
  indexName,
  txpData,
  source,
  organism,
  release,
  genome
)
```

## Arguments

- digest:

  character string of the full digest of the reference transcripts, see
  [`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
  with `fullDigest=TRUE`

- digestType:

  character string of the digest, default `"sha256"`

- indexName:

  a name for the `index` when storing the linkedTxpData,

- txpData:

  *GRanges* providing information about ranges representing the
  transcript sequences linked to `digest`

- source:

  the source of transcriptome, e.g. `denovo`. See
  [`makeLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md)
  for more information on specifying source

- organism:

  organism (e.g. "Homo sapiens")

- release:

  release number (e.g. "27")

- genome:

  genome (e.g. "GRCh38", or "none")

## Value

nothing, the function is run for its side effects

## Details

The `txpData` object is saved in the
[`getTximetaBFC()`](https://thelovelab.github.io/tximeta/reference/getTximetaBFC.md)
location, appended with a 32-character substring of `digest`. The tibble
listing all `linkedTxpData` is named `linkedTxpDataTbl` and is listed in
the same location.

## Examples

``` r

novel <- data.frame(seqnames = paste0("chr", rep(1:22, each=500)),
  start = 1e6 + 1 + 0:499 * 1000, end = 1e6 + 1 + 0:499 * 1000 + 1000 - 1,
  strand = "+", tx_name = paste0("novel", 1:(22*500)), 
  gene_id = paste0("novel_gene", rep(1:(22*10), each=50)), 
type = "protein_coding")
novel_gr <- as(novel, "GRanges")
names(novel_gr) <- novel$tx_name

makeLinkedTxpData(
 digest = "43158f2c8e88e3acd77c22aee557625a6f1b6a5038cfc7deb5e64903892d8070",
 digestType = "sha256",
 indexName = "my_novel_txps",
 txpData = novel_gr,
 source = "novel", organism="Homo sapiens", 
 release="v1", genome="GRCh38"
)
#> linking user-provided metadata to digest: 43158f...
#> saving txpData object in bfc
#> saving linkedTxpData in bfc (first time)

# to clear the entire linkedTxome table
# (don't run unless you want to clear this table!)
# bfcloc <- getTximetaBFC()
# bfc <- BiocFileCache(bfcloc)
# bfcremove(bfc, bfcquery(bfc, "linkedTxpDataTbl")$rid)
```
