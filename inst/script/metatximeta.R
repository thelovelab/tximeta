hashit <- function(source, organism, release) {
  if (source == "Gencode") {
    if (organism == "Homo sapiens") {
      fasta <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                      release,"/gencode.v",release,".transcripts.fa.gz")
      gtf <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
                    release,"/gencode.v",release,".annotation.gtf.gz")
    }
    if (organism == "Mus musculus") {
      fasta <- paste0("ftp://ftp.sanger.ac.uk/pub/databases/gencode/Gencode_mouse/release_",
                      release,"/gencode.v",release,".transcripts.fa.gz")
      gtf <- paste0("ftp://ftp.sanger.ac.uk/pub/databases/gencode/Gencode_mouse/release_",
                    release,"/gencode.v",release,".annotation.gtf.gz")
    }
  }
  if (source == "Ensembl") {
    org <- sub(" ","_",tolower(organism))
    org_upper <- sub(" ","_",organism)
    fasta <- paste0("ftp://ftp.ensembl.org/pub/release-",
                    release,"/fasta/",org,"/cdna/",org_upper,".GRCh38.cdna.all.fa.gz")
    gtf <- paste0("ftp://ftp.ensembl.org/pub/release-",
                  release,"/gtf/",org,"/",org_upper,".GRCh38.",release,".gtf.gz")
  }
  download.file(fasta, "transcript.fa")
  system("compute_fasta_digest")
}
