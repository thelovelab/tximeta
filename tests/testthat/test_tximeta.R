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

  expect_warning({se <- tximeta(coldata)}, "missing from the GTF")
  
})

test_that("tximeta can import kallisto", {

  library(readr)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"kallisto", samples$run, "abundance.tsv.gz")
  coldata <- data.frame(files, names=paste0("sample",1:6))
  se <- tximeta(coldata, type="kallisto", txOut=TRUE)
  
})
