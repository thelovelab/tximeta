# catNC = whether to concatenate the non-coding RNA to the coding (Ensembl only)
hashit <- function(source, organism, release, catNC=FALSE) {

  if (source == "Gencode") {
    if (organism == "Homo sapiens") {
      fasta <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                      release,"/gencode.v",release,".transcripts.fa.gz")
      gtf <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                    release,"/gencode.v",release,".annotation.gtf.gz")
      genome <- "GRCh38"
    }
    if (organism == "Mus musculus") {
      fasta <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_",
                      release,"/gencode.v",release,".transcripts.fa.gz")
      gtf <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_",
                    release,"/gencode.v",release,".annotation.gtf.gz")
      genome <- "GRCm38"
    }
  }
  
  if (source == "Ensembl") {
    org <- sub(" ","_",tolower(organism))
    org_upper <- sub(" ","_",organism)
    genome <- if (organism == "Homo sapiens") {
                "GRCh38"
              } else if (organism == "Mus musculus") {
                  "GRCm38"
              } else if (organism == "Drosophila melanogaster") "BDGP6"
    fasta <- paste0("ftp://ftp.ensembl.org/pub/release-",
                    release,"/fasta/",org,"/cdna/",org_upper,".",genome,".cdna.all.fa.gz")
    if (catNC) {
      nc <- paste0("ftp://ftp.ensembl.org/pub/release-",
                   release,"/fasta/",org,"/ncrna/",org_upper,".",genome,".ncrna.fa.gz")
      fasta <- c(fasta, nc)
    }
    gtf <- paste0("ftp://ftp.ensembl.org/pub/release-",
                  release,"/gtf/",org,"/",org_upper,".",genome,".",release,".gtf.gz")

  }

  if (source == "RefSeq") {
    genome <- if (organism == "Homo sapiens") {
                "GRCh38"
              } else if (organism == "Mus musculus") {
                "GRCm38"
              }
    fasta <- ""
    gtf <- ""
  }

  stopifnot(all(sapply(fasta, RCurl::url.exists)))
  stopifnot(RCurl::url.exists(gtf))
  if (catNC) {
    for (f in fasta) {
      download.file(f, basename(f))
    }
    system(paste("cat",paste(basename(fasta),collapse=" "),"> transcripts.fa.gz"))
    for (f in fasta) {
      system(paste("rm -f",basename(f)))
    }
  } else {
    download.file(fasta, "transcripts.fa.gz")
  }
  system("gunzip transcripts.fa.gz")
  system("compute_fasta_digest --reference transcripts.fa --out hash.json")
  sha256 <- jsonlite::fromJSON("hash.json")$seq_hash
  df <- data.frame(source, organism, release, genome, paste(fasta, collapse=" "), gtf, sha256)
  write.table(df, file="newrows.csv", append=TRUE, quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)
  system("rm -f transcripts.fa")
  system("rm -f hash.json")
  NULL
}

hashit("Gencode", "Homo sapiens", "30")
hashit("Gencode", "Mus musculus", "M21")

hashit("Ensembl", "Homo sapiens", "95")
hashit("Ensembl", "Homo sapiens", "96")
hashit("Ensembl", "Mus musculus", "95")
hashit("Ensembl", "Mus musculus", "96")
hashit("Ensembl", "Drosophila melanogaster", "95")
# hashit("Ensembl", "Drosophila melanogaster", "96") => genome moved from BDGP6 to BDGP6.22

hashit("Ensembl", "Homo sapiens", "95", catNC=TRUE)
hashit("Ensembl", "Homo sapiens", "96", catNC=TRUE)
hashit("Ensembl", "Mus musculus", "95", catNC=TRUE)
hashit("Ensembl", "Mus musculus", "96", catNC=TRUE)

