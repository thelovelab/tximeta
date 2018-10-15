hashit <- function(source, organism, release) {
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
                "GRCh38" }
              else if (organism == "Mus musculus") {
                  "GRCm38" }
              else if (organism == "Drosophila melanogaster") "BDGP6"
    fasta <- paste0("ftp://ftp.ensembl.org/pub/release-",
                    release,"/fasta/",org,"/cdna/",org_upper,".",genome,".cdna.all.fa.gz")
    gtf <- paste0("ftp://ftp.ensembl.org/pub/release-",
                  release,"/gtf/",org,"/",org_upper,".",genome,".",release,".gtf.gz")

  }
  stopifnot(RCurl::url.exists(fasta))
  stopifnot(RCurl::url.exists(gtf))
  download.file(fasta, "transcripts.fa.gz")
  system("gunzip transcripts.fa.gz")
  system("compute_fasta_digest --reference transcripts.fa --out hash.json")
  sha256 <- jsonlite::fromJSON("hash.json")$seq_hash
  df <- data.frame(source, organism, release, genome, fasta, gtf, sha256)
  write.table(df, file="newrows.csv", append=TRUE, quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)
  system("rm -f transcripts.fa")
  system("rm -f hash.json")
  NULL
}

hashit("Gencode", "Homo sapiens", "29")
hashit("Gencode", "Mus musculus", "M19")
hashit("Gencode", "Mus musculus", "M18")
hashit("Ensembl", "Homo sapiens", "94")
hashit("Ensembl", "Homo sapiens", "93")
hashit("Ensembl", "Mus musculus", "94")
hashit("Ensembl", "Mus musculus", "93")
hashit("Ensembl", "Drosophila melanogaster", "94")
hashit("Ensembl", "Drosophila melanogaster", "93")
