# Creation of the Salmon indices in `extdata`

The `extdata` directory contains two directories, each with a
`header.json` file:

* `Drosophila_melanogaster.BDGP6.cdna.v92_salmon_0.10.2`
* `Drosophila_melanogaster.BDGP6.v92_salmon_0.10.2`

These directories were created by running the `index` step of *Salmon*
version 0.10.2 on the Drosophila melanogaster transcriptome. The only
files that have been posted are the `header.json` files from the index
directory, as the index itself is not needed by `tximeta` and is too
large to include. These directories allow us to demonstrate the use of
`tximeta` when the transcriptome is not known.

The transcriptome FASTA files were downloaded from the following FTP sites:

* <ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz>
* <ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.ncrna.fa.gz>

The `Drosophila_melanogaster.BDGP6.v92_salmon_0.10.2` index was
created by `cat`-ing the `cdna` and `ncrna` gzipped FASTA together
before running `salmon index`.
