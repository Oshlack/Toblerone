# Toblerone

Toblerone is a method for detecting internal exon (not first or last) deletions in RNA-seq data. It is currently in development, extending the methods described here: https://github.com/Oshlack/ALL-RNAseq-utility-paper


## Strucuture


### Pseudoaligner

The psuedaligner core is currently called `tinyt` and can be run in standalone mode form the command line, or is inclued in the Toblerone app (see below).

```
tinyt

Usage:
  tinyt index [--num-threads=<n>] -i <index> <ref-fasta>
  tinyt map [--num-threads=<n>] [--read-length=<r>][--trim-size=<t>] [--skip-trim] [--output=<file>] -i <index> <reads-fastq> [<reads-pair-fastq>]
  tinyt -h | --help | -v | --version

Options:
  -n --num-threads N  Number of worker threads [default: 2]
  -t --trim-size T    Size of base pairs to trim when checking unique read matches [default: 2]
  -s --skip-trim      Skip the trim read check for unqiue read matches
  -r --read-length R  Provide read length for depth estimation
  -o --output FILE    Output results to file instead of stdout
  -h --help           Show this screen.
  -v --version        Show version.
```

#### Index

For a candidate gene, e.g. IKZF1, we take the the canonical transcript and generate a specialized transcriptome reference that consists of the original transcript plus deletion transcripts. Deletion transcripts consist of combinations of continuous exon deletion, excluding edge exons (first and last) exons. For N exons, (N-1 choose 2) additional deletion transcripts created.

#####  Generate transcripts


##### Create index



#### Map


#### Code & License

This pseudoalignment implementation is based on the 10X Genomics Pseudoaligner code (https://github.com/10XGenomics/rust-pseudoaligner/), which itself draws on the concepts from Kalisto(), Salmon() and others. It is released under the MIT license in line with the template source. It is heavily modified for a Toblerone index, and not intended for general transcriptomes.  



### App

TBC


# Roadmap

- [ ] Move transcriptome generation script into core 

# Acknowledgements 


