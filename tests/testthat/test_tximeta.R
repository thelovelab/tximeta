context("tximeta")
test_that("tximeta works as expected", {

  dir <- system.file("extdata/salmon_dm", package="tximportData")
  files <- file.path(dir, "SRR1197474_cdna", "quant.sf.gz") 
  coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
  
  dir <- system.file("extdata", package="tximeta")
  indexDir <- file.path(dir, "Drosophila_melanogaster.BDGP6.cdna.v92_salmon_0.10.2")
  fastaFTP <- "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
  
  dir2 <- system.file("extdata/salmon_dm", package="tximportData")
  gtfPath <- file.path(dir2,"Drosophila_melanogaster.BDGP6.92.gtf.gz")
  makeLinkedTxome(indexDir=indexDir, source="Ensembl", organism="Drosophila melanogaster",
                  release="92", genome="BDGP6", fasta=fastaFTP, gtf=gtfPath, write=FALSE)

  expect_warning({se <- tximeta(coldata)}, "annotation is missing")

  # check adding exons
  library(SummarizedExperiment)
  se <- addExons(se)

  # check summarize to gene
  gse <- summarizeToGene(se)
  expect_error(addExons(gse), "transcript-level")

  # check adding IDs
  library(org.Dm.eg.db)
  gse <- addIds(gse, "REFSEQ", gene=TRUE)
    
  # just a vector of file paths is ok
  expect_warning({se <- tximeta(files)})

  # check error on txOut=FALSE
  expect_error({se <- tximeta(coldata, txOut=FALSE)}, "transcript-level output")

  # check skipping metadata, appropriate warnings
  se <- tximeta(coldata, skipMeta=TRUE)
  expect_error({summarizeToGene(se)}, "transcriptome metadata")
  expect_error({addIds(se)}, "transcriptome metadata")
  
})

test_that("tximeta can import GENCODE and Ensembl", {

  if (FALSE) {
    dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
    coldata <- data.frame(files, names=paste0("sample",1:6))
    se <- tximeta(coldata)

    dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    files <- file.path(dir,"salmon_gibbs", samples$run, "quant.sf.gz")
    coldata <- data.frame(files, names=paste0("sample",1:6))
    se <- tximeta(coldata, dropInfReps=TRUE)
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
  
})
