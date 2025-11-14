# Package index

## Package overview

- [`tximeta-package`](https://thelovelab.github.io/tximeta/reference/tximeta-package.md)
  : Import transcript quantification with metadata

## Main import functions

- [`tximeta()`](https://thelovelab.github.io/tximeta/reference/tximeta.md)
  : Import transcript quantification with metadata
- [`summarizeToGene(`*`<SummarizedExperiment>`*`)`](https://thelovelab.github.io/tximeta/reference/summarizeToGene.md)
  : Summarize estimated quantitites to gene-level

## linkedTxome: linking data to metadata

- [`makeLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md)
  [`loadLinkedTxome()`](https://thelovelab.github.io/tximeta/reference/linkedTxome.md)
  : Make and load linked transcriptomes (linked GTF and FASTA)
- [`makeLinkedTxpData()`](https://thelovelab.github.io/tximeta/reference/linkedTxpData.md)
  : Make linked transcript data (linked GRanges)

## Mixed reference transcripts

- [`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md)
  : Import quantification across mixed reference transcripts

- [`inspectDigests()`](https://thelovelab.github.io/tximeta/reference/inspectDigests.md)
  :

  Inspect digest matches from
  [`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md)
  imported data

- [`updateMetadata()`](https://thelovelab.github.io/tximeta/reference/updateMetadata.md)
  :

  Update transcript metadatda for
  [`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md)
  imported data

## Add information to object

- [`addCDS()`](https://thelovelab.github.io/tximeta/reference/addCDS.md)
  : Add CDS to rowRanges of a transcript-level SummarizedExperiment
- [`addExons()`](https://thelovelab.github.io/tximeta/reference/addExons.md)
  : Add exons to rowRanges of a transcript-level SummarizedExperiment
- [`addIds()`](https://thelovelab.github.io/tximeta/reference/addIds.md)
  : Add IDs to rowRanges of a SummarizedExperiment

## Retrieve information from object

- [`retrieveCDNA()`](https://thelovelab.github.io/tximeta/reference/retrieveCDNA.md)
  : Retrieve the cDNA transcript sequence for a SummarizedExperiment
- [`retrieveDb()`](https://thelovelab.github.io/tximeta/reference/retrieveDb.md)
  : Retrieve the TxDb or EnsDb associated with a SummarizedExperiment

## Other helpers

- [`splitSE()`](https://thelovelab.github.io/tximeta/reference/splitSE.md)
  : Split SummarizedExperiment by gene categories
- [`makeDGEList()`](https://thelovelab.github.io/tximeta/reference/makeDGEList.md)
  : Make a DGEList from tximeta output
- [`getTximetaBFC()`](https://thelovelab.github.io/tximeta/reference/getTximetaBFC.md)
  [`setTximetaBFC()`](https://thelovelab.github.io/tximeta/reference/getTximetaBFC.md)
  : Get or set the directory of the BiocFileCache used by tximeta
