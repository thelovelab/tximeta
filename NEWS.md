# tximeta 1.23.1

* GENCODE 46 (H.s.), M35 (M.m), and Ensembl 112

# tximeta 1.21.4

* Changing language in docs to "digest" instead of "checksum".

# tximeta 1.21.3

* GENCODE 44 (H.s.), M34 (M.m), and Ensembl 111
* RefSeq p13 for human, p6 for mouse

# tximeta 1.20.0

* Add argument to summarizeToGene(): `assignRanges` that takes
  either "range" (default) or "abundant", and determines
  the ranges that are attached to the SE (rowRanges).
  Note that this new argument does not affect the data aggregation
  at all (counts and abundance are summarized to gene by `tximport`).
  The default behavior of summarizeToGene() returns ranges that
  correspond to the `range` of the isoforms of the gene, that is the
  leftmost basepair to the rightmost basepair of any isoform.
  The non-default "abundant" instead returns the range of the
  most abundant isoform in the data, averaging over samples.
  Information about the choice of range is included in `mcols`
* Added support for piscem-infer: `type="piscem"` also 
  auto-detected from file ending.

# tximeta 1.19.8

* Added support for piscem-infer: `type="piscem"` also 
  auto-detected from file ending.

# tximeta 1.19.6

* Fixed genome build for mouse M26 and higher to GRCm39,
  thanks to Charlotte Soneson.

# tximeta 1.14.0

* Allow GTF specification in linkedTxome to be a serialized
  GRanges file (a file path to a .rda or .RData file). This
  bypasses some apparent issue where makeTxDbFromGFF fails
  while makeTxDbFromGRanges works.
* Up to GENCODE 40 (H.s.), M29 (M.m), and Ensembl 106

# tximeta 1.10.0

* Added more tximeta() messaging about specifying the
  'source' in linkedTxome. Essentially, this triggers
  GTF processing behavior that users may want to avoid,
  and so specifying a string other than "Ensembl" may be
  preferred. Also added note to vignette.
* Added note to vignette about alevin import with tximeta
  where the 'tgMap' step requires gene IDs and not gene
  symbols.
* Added hashes for:
  GENCODE 38 (H.s.), M27 (M.m), and Ensembl 104;
  GENCODE 37 (H.s.), M26 (M.m), and Ensembl 103;
  GENCODE 36 and Ensembl 102.
* Fixed a bug around AnnotationHub pulldown when using RefSeq
  as the source.
* Fixed a bug where multiple parsed Ensembl GTF TxDb would be
  added to the BiocFileCache with the same rname.

# tximeta 1.8.0

* Added 'fromDb' argument to addIds() to allow IDs to be
  added from the associated TxDb/EnsDb instead of the org
  package (which is used by default).
  Feature suggestion from Kristoffer Vitting-Seerup.
* Added function retrieveCDNA() that will download or load
  a cached version of the transcript sequences used for
  quantification. Note that the returned sequences are not
  ordered or matched to the rows of the SummarizedExperiment
  object. Feature suggestion from Kristoffer Vitting-Seerup.
* Added function addCDS() that will add CDS ranges for coding
  transcripts (and fills in original ranges for non-coding),
  as well as a 'coding' column as a logical indicator.
  Feature suggestion from Kristoffer Vitting-Seerup.
* Added option that environmental variable TXIMETA_HUB_CACHE
  can be used to set tximeta's cache location, to avoid
  prompting the user on the first run of tximeta().
* tximeta() will now pull GENCODE TxDb from AnnotationHub
  when it is listed there (only Homo sapiens are at this
  point in time). Thanks to Leonardo Collado-Torres for
  the suggestion!
* Now summarizeToGene() will add a column tx_ids, which is a
  CharacterList of the transcript IDs.

# tximeta 1.6.0

* Added PLOS Computational Biology citation! :-)
* Added function `splitSE` to split one assay of a 
  SummarizedExperiment into multiple assays, each containing 
  features of a given type.
* Added a wrapper function makeDGEList() to simplify making
  a DGEList for use with edgeR. See vignette for example.
* tximeta will now make use of EnsDb created and distributed
  on AnnotationHub, unless useHub=FALSE. Also, a new function
  retrieveDb() can be called on a SummarizedExperiment to
  retrieve the underlying TxDb or EnsDb.
* tximeta can now use `customMetaInfo` argument to locate
  a custom metadata information file such as `meta_info.json`,
  which should contain a tage, `index_seq_hash`, with the SHA-256
  hash value of the reference transcripts.
* Added `markDuplicateTxps` argument to add `hasDuplicate`
  and `duplicates` columns to rowData of SummarizedExperiment.
  One note is that, for efficiency, this argument and
  `cleanDuplicateTxps` will now share a duplicates CharacterList
  that is stored in the BiocFileCache, with the name `dups-...`.
  Therefore, if you have previously used `cleanDuplicateTxps`,
  you may need to bfcremove() any `dups-...` entries.
  Summarization to gene level will keep track of `numDupSets`
  per gene which informs about the number of transcripts sets
  (equivalence classes by transcript sequence content).
* If during the indexing step, user didn't use --gencode
  for a Gencode transcriptome file, tximeta will deal with
  this internally now by stripping all characters after the
  vertical bar `|`, in order to match long transcript names
  in the `quant.sf` files to the correct transcript names
  in the GTF.

# tximeta 1.4.0

* tximeta will now pull down RefSeq seqinfo, using the
  dirname() of the GTF location, and assuming some
  consistency in the structure of the assembly_report.txt
  that is located in the same directory. Needs more
  testing though across releases and organisms.
* expanded caching of ranges to exons and genes as well.
  Exons in particular take a long time to build from
  TxDb, so this saves quite a lot of time.
* new 'addExons' function will add exons to trancript-level
  summarized experiments, by replacing transcript GRanges
  with exon-by-transcript GRangesList. Purposely designed
  only for transcript-level, see note in ?addExons
* tximeta now also caches the transcript ranges themselves,
  rather than just the TxDb. This shaves extra seconds off
  the tximeta() call!
* add 'skipSeqinfo' argument, which avoids attempting
  to fetch chromosome information (from UCSC) if set
  to TRUE.

# tximeta 1.2.0

* Specifying gene=TRUE in addIds() when rows are
  transcripts will attempt to use a gene_id column
  to map the IDs. This usually gives a better mapping
  rate.
* Cut off version number from Ensembl names only (not GENCODE)
* Added 'cleanDuplicateTxps' argument, which does a lot of work
  for the user: it downloads the FASTA from the source, identifies
  duplicate transcripts (identical cDNA sequence) then looks to
  see if transcripts that are in the quantification files,
  but missing from the GTF, could be renamed from the list of
  duplicate transcripts such that they would be present in the GTF.
* Added coding * non-coding combinations of Ensembl transcriptomes
  to the hash table. Must be in this order: coding, then non-coding.
* Added support for dammit de novo transcriptomes.
* Added summarizeToGene as a method, to avoid conflicts with tximport.
* Added in Charlotte's code to split out GENCODE and Ensembl code
  for generating transcript ranges.
