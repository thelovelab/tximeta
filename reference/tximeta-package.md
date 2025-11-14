# Import transcript quantification with metadata

The tximeta package imports abundances (TPM), estimated counts, and
effective lengths from quantification tools, and will output a
*SummarizedExperiment* (SE) object. For salmon and related
quantification tools,
[`tximeta()`](https://thelovelab.github.io/tximeta/reference/tximeta.md)
will attempt to identify the correct provenance of the reference
transcripts and automatically attach the transcript ranges to the
SummarizedExperiment, to facilitate downstream integration with other
datasets. The automatic identification of reference transcripts should
work out-of-the-box for human or mouse transcriptomes from the sources:
GENCODE, Ensembl, or RefSeq. See also
[`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md)
for importing data when the reference transcripts were derived from a
mix of annotated (e.g. GENCODE) and novel or custom transcripts.

## Details

The main functions are:

- [`tximeta()`](https://thelovelab.github.io/tximeta/reference/tximeta.md) -
  with key argument `coldata` specifying sample information

- [`summarizeToGene()`](https://thelovelab.github.io/tximeta/reference/summarizeToGene.md) -
  summarize quantification to gene-level

- [`importData()`](https://thelovelab.github.io/tximeta/reference/importData.md) -
  import quantification with mixed reference transcript sets

All software-related questions should be posted to the Bioconductor
Support Site:

<https://support.bioconductor.org>

The code can be viewed at the GitHub repository, which also lists the
contributor code of conduct:

<https://github.com/thelovelab/tximeta>

## References

- *tximeta* reference:

Michael I. Love, Charlotte Soneson, Peter F. Hickey, Lisa K. Johnson N.
Tessa Pierce, Lori Shepherd, Martin Morgan, Rob Patro (2020) *Tximeta:
reference sequence checksums for provenance identification in RNA-seq*.
PLOS Computational Biology.
<https://doi.org/10.1371/journal.pcbi.1007664>

- *tximport* reference (the effective length GLM offset and
  counts-from-abundance):

Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015)
*Differential analyses for RNA-seq: transcript-level estimates improve
gene-level inferences*. F1000Research.
<http://doi.org/10.12688/f1000research.7563>

## See also

Useful links:

- <https://github.com/thelovelab/tximeta>

## Author

Michael I. Love, Charlotte Soneson, Peter Hickey, Rob Patro
