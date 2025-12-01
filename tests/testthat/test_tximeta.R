context("tximeta")
test_that("tximeta works as expected", {

  dir <- system.file("extdata/salmon_dm", package="tximportData")
  files <- file.path(dir, "SRR1197474", "quant.sf") 
  coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
  
  indexDir <- file.path(dir, "Dm.BDGP6.22.98_salmon-0.14.1")
  fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
                "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")
  
  gtfPath <- file.path(dir,"Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
  makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster",
                  release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)

  # TODO not throwing warnings on Bioc machines?
  expect_warning({se <- tximeta(coldata)}, "the annotation is missing")
  #se <- tximeta(coldata)

  # check adding IDs
  library(org.Dm.eg.db)
  se <- addIds(se, column="SYMBOL")
  
  # check adding IDs from TxDb/EnsDb (commented out due to Ahub/EnsDb issue)
  # se <- addIds(se, column="SYMBOL", fromDb=TRUE)
  
  # check adding CDS
  library(SummarizedExperiment)
  se_cds <- addCDS(se)
  
  # check adding exons
  se_exons <- addExons(se)

  # check adding CDS on exons
  se_cds <- addCDS(se)

  # check summarize to gene
  gse <- summarizeToGene(se)
  expect_error(addExons(gse), "transcript-level")

  gseAbndt <- summarizeToGene(se, assignRanges="abundant")

  gseAbndt <- summarizeToGene(se_exons, assignRanges="abundant")
  
  # check adding IDs
  library(org.Dm.eg.db)
  gse <- addIds(gse, "REFSEQ", gene=TRUE)

  # check retrieving the database
  edb <- retrieveDb(se)

  # check retrieving the cDNA sequence
  if (FALSE) {
    # this requires access to Ensembl ftp
    cdna <- retrieveCDNA(se)
  }
  
  expect_warning({se <- tximeta(files)}, "the annotation is missing")

  # check error on txOut=FALSE
  expect_error({se <- tximeta(coldata, txOut=FALSE)},
               "transcript-level output")

  # check skipping metadata, appropriate warnings
  se <- tximeta(coldata, skipMeta=TRUE)
  expect_error({summarizeToGene(se)}, "transcriptome metadata")
  expect_error({addIds(se)}, "transcriptome metadata")

  # check customMetaInfo
  expect_warning(
    {se <- tximeta(coldata, customMetaInfo="aux_info/meta_info.json")}, 
    "the annotation is missing"
  )
  expect_error(tximeta(coldata, customMetaInfo="foobar.json"),
               "metadata files are missing")

  # check make DGEList
  library(edgeR)
  y <- makeDGEList(se)
  
})

test_that("tximeta can import GENCODE and Ensembl", {

  # breaks with no internet
  if (FALSE) {

    ### GENCODE ###
    dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
    coldata <- data.frame(files, names=paste0("sample",1:6))

    # with AnnotationHub (default)
    se <- tximeta(coldata)
    expect_true(metadata(se)$txomeInfo$source == "GENCODE")
    
    gse <- summarizeToGene(se)

    # check adding IDs from TxDb/EnsDb
    se <- addIds(se, column="TXTYPE", fromDb=TRUE)
    
    # without AnnotationHub
    se <- tximeta(coldata, useHub=FALSE)

    ### Ensembl ###
    dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    files <- file.path(dir,"salmon_gibbs", samples$run, "quant.sf.gz")
    coldata <- data.frame(files, names=paste0("sample",1:6))

    # with AnnotationHub (default)
    se <- tximeta(coldata, dropInfReps=TRUE, useHub=TRUE)
    
    # without AnnotationHub
    se <- tximeta(coldata, dropInfReps=TRUE, useHub=FALSE)

    gse <- summarizeToGene(se)

    # test the mark duplicate code
    
    dir <- system.file("extdata", package="macrophage")
    samples <- list.files(file.path(dir, "quants"))[1:3]
    files <- file.path(dir,"quants", samples, "quant.sf.gz")
    coldata <- data.frame(files, names=paste0("sample",1:3))

    se <- tximeta(coldata, dropInfReps=TRUE, markDuplicateTxps=TRUE)
    gse <- summarizeToGene(se)
    
  }

})

test_that("tximeta can import inferential replicates", {

  # breaks with no internet
  if (FALSE) {
    library(SummarizedExperiment)

    # check the GEUVADIS samples with Salmon Gibbs samples
    dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    files <- file.path(dir,"salmon_gibbs", samples$run, "quant.sf.gz")
    coldata <- data.frame(files, names=paste0("sample",1:6))
    se <- tximeta(coldata)
    expect_true("infRep1" %in% assayNames(se))

    se2 <- tximeta(coldata, varReduce=TRUE)
    expect_true("variance" %in% assayNames(se2))

    gse <- summarizeToGene(se)
    expect_true("infRep1" %in% assayNames(gse))

    gse <- summarizeToGene(se, varReduce=TRUE)
    expect_true("variance" %in% assayNames(gse))

    # check the macrophage dataset with Salmon Gibbs samples

    dir <- system.file("extdata", package="macrophage")
    coldata <- read.csv(file.path(dir, "coldata.csv"))
    coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
    coldata <- coldata[1:2,]
    se <- tximeta(coldata)
    se <- tximeta(coldata, skipSeqinfo=TRUE)
    
  }
  
})

# refseq GTF files can't be imported by txdbmaker::makeTxDbFromGFF (August 2025)...
# so then tximeta won't really work anymore

# test_that("tximeta can import refseq", {
#     dir <- system.file("extdata", package="tximportData")
#     files <- file.path(dir,"refseq/ERR188021/quant.sf.gz")
#     file.exists(files)
#     coldata <- data.frame(files, names="A")
#     # here the important test is if we can pull down seqinfo
#     # (sequence names and lengths) from RefSeq FTP site
#     se <- tximeta(coldata)
#     seqinfo(se)
# })

test_that("tximeta can import kallisto", {

  # test requires rhdf5...
  if (FALSE) {
    dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    files <- file.path(dir,"kallisto", samples$run, "abundance.tsv.gz")
    coldata <- data.frame(files, names=paste0("sample",1:6))
    se <- tximeta(coldata, type="kallisto", txOut=TRUE)
    
    # inferential replicates as well
    library(SummarizedExperiment)
    files <- file.path(dir,"kallisto_boot", samples$run, "abundance.h5")
    coldata <- data.frame(files, names=paste0("sample",1:6))
    se <- tximeta(coldata, type="kallisto", txOut=TRUE)
    expect_true("infRep1" %in% assayNames(se))
    se <- tximeta(coldata, type="kallisto", txOut=TRUE, varReduce=TRUE)
    expect_true("variance" %in% assayNames(se))
  }
  
})
