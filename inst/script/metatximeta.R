hashit <- function(source, organism, release, catNC=FALSE, ...) {

  if (source == "GENCODE") {
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
      if (as.numeric(sub("M", "", release)) >= 26) {
        genome <- "GRCm39"
      } else {
        genome <- "GRCm38"
      }
    }
  }
  
  if (source == "Ensembl") {
    org <- sub(" ","_",tolower(organism))
    org_upper <- sub(" ","_",organism)
    genome <- if (organism == "Homo sapiens") {
                "GRCh38"
              } else if (organism == "Mus musculus") {
                if (as.numeric(release) >= 103) {
                  "GRCm39"
                } else {
                  "GRCm38"
                }
              } else if (organism == "Drosophila melanogaster") {
                if (as.numeric(release) >= 110) {
                  "BDGP6.46"
                } else if (as.numeric(release) >= 103) {
                  "BDGP6.32"
                } else if (as.numeric(release) >= 99) {
                  "BDGP6.28"
                } else if (as.numeric(release) >= 96) {
                  "BDGP6.22"
                } else {
                  "BDGP6"
                }
              } 
  
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
    org_upper <- sub(" ","_",organism)
    chop <- function(x) as.numeric(sub(".*\\.p","",x))
    if (organism == "Homo sapiens") {
      genome <- "GRCh38"
      release_num <- chop(release) + 26
      gcf <- "GCF_000001405"
    } else if (organism == "Mus musculus") {
      genome <- "GRCm38"
      release_num <- chop(release) + 20
      gcf <- "GCF_000001635"
    }
    fasta <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/",org_upper,
                    "/all_assembly_versions/",gcf,".",release_num,"_",release,"/",
                    gcf,".",release_num,"_",release,"_rna.fna.gz")
    gtf <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/",org_upper,
                  "/all_assembly_versions/",gcf,".",release_num,"_",release,"/",
                  gcf,".",release_num,"_",release,"_genomic.gff.gz")
  }

  #message(fasta)
  #message(gtf)
  
  stopifnot(all(sapply(fasta, RCurl::url.exists)))
  stopifnot(RCurl::url.exists(gtf))
  if (catNC) {
    for (f in fasta) {
      download.file(f, basename(f), ...)
    }
    system(paste("cat",paste(basename(fasta),collapse=" "),"> transcripts.fa.gz"))
    for (f in fasta) {
      system(paste("rm -f",basename(f)))
    }
  } else {
    download.file(fasta, "transcripts.fa.gz", ...)
  }
  system("gunzip transcripts.fa.gz")
  system("compute_fasta_digest --reference transcripts.fa --out hash.json")
  sha256 <- jsonlite::fromJSON("hash.json")$seq_hash
  df <- data.frame(source, organism, release, genome, paste(fasta, collapse=" "), gtf, sha256)
  write.table(df, file="newrows.csv", append=TRUE, quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)
  system("rm -f transcripts.fa")
  system("rm -f hash.json")
  message("done!")
}

# download method
m <- "wget"

i <- 45
hashit("GENCODE", "Homo sapiens", i, method=m)
i <- 34
hashit("GENCODE", "Mus musculus", paste0("M",i), method=m)
#
i <- 111
hashit("Ensembl", "Homo sapiens", i, method=m)
hashit("Ensembl", "Mus musculus", i, method=m)
hashit("Ensembl", "Drosophila melanogaster", i, method=m)
#
hashit("Ensembl", "Homo sapiens", i, catNC=TRUE, method=m)
hashit("Ensembl", "Mus musculus", i, catNC=TRUE, method=m)
hashit("Ensembl", "Drosophila melanogaster", i, catNC=TRUE, method=m)
#
# RefSeq permalink transcriptomes are rarely updated
hashit("RefSeq", "Homo sapiens", paste0("GRCh38.p",i), method=m)
hashit("RefSeq", "Mus musculus", paste0("GRCm38.p",i), method=m)
