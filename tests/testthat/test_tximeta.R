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

  # TODO why not throwing warnings on Bioc
  #expect_warning({se <- tximeta(coldata)}, "annotation is missing")
  se <- tximeta(coldata)

  # check adding IDs
  library(org.Dm.eg.db)
  se <- addIds(se, column="SYMBOL")
  
  # check adding IDs from TxDb/EnsDb (commented out due to Ahub/EnsDb issue)
  # se <- addIds(se, column="SYMBOL", fromDb=TRUE)
  
  # check adding CDS
  library(SummarizedExperiment)
  se2 <- addCDS(se)
  
  # check adding exons
  se <- addExons(se)

  # check adding CDS on exons
  se <- addCDS(se)

  # check summarize to gene
  gse <- summarizeToGene(se)
  expect_error(addExons(gse), "transcript-level")

  # check adding IDs
  library(org.Dm.eg.db)
  gse <- addIds(gse, "REFSEQ", gene=TRUE)

  # check retrieving the database
  edb <- retrieveDb(se)

  # check retrieving the cDNA sequence
  cdna <- retrieveCDNA(se)
  
  # just a vector of file paths is ok...
  # TODO why not throwing warnings on Bioc
  #expect_warning({se <- tximeta(files)})
  se <- tximeta(files)

  # check error on txOut=FALSE
  expect_error({se <- tximeta(coldata, txOut=FALSE)},
               "transcript-level output")

  # check skipping metadata, appropriate warnings
  se <- tximeta(coldata, skipMeta=TRUE)
  expect_error({summarizeToGene(se)}, "transcriptome metadata")
  expect_error({addIds(se)}, "transcriptome metadata")

  # check customMetaInfo
  se <- tximeta(coldata, customMetaInfo="aux_info/meta_info.json")
  expect_error(tximeta(coldata, customMetaInfo="foobar.json"),
               "metadata files are missing")

  # check make DGEList
  library(edgeR)
  y <- makeDGEList(se)
  
})

test_that("tximeta can import GENCODE and Ensembl", {

  if (FALSE) {

    ### GENCODE ###
    dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
    coldata <- data.frame(files, names=paste0("sample",1:6))

    # with AnnotationHub (default)
    se <- tximeta(coldata)
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

  # don't want to rely on internet connection for tests...
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

test_that("tximeta can import alevin", {

  dir <- system.file("extdata", package="tximportData")
  files <- file.path(dir,"alevin/neurons_900_v014/alevin/quants_mat.gz")
  file.exists(files)
  coldata <- data.frame(files, names="neurons")
  #se <- tximeta(coldata, type="alevin")
  #se <- tximeta(coldata, type="alevin", dropInfReps=TRUE)
  #se <- tximeta(coldata, type="alevin", dropInfReps=TRUE, skipMeta=TRUE)
  
})

test_that("tximeta can import refseq", {

  if (FALSE) {
    dir <- system.file("extdata", package="tximportData")
    files <- file.path(dir,"refseq/ERR188021/quant.sf.gz")
    file.exists(files)
    coldata <- data.frame(files, names="A")

    # here the important test is if we can pull down seqinfo
    # (sequence names and lengths) from RefSeq FTP site
    se <- tximeta(coldata)
    seqinfo(se)
  }
    
})

test_that("tximeta can import kallisto", {

  if (requireNamespace(package="rhdf5", quietly=TRUE)) {
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
