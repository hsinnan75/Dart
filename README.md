DART: a fast and accurate RNA-seq mapper with a divide and conquer strategy
===================

Developers: Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu Institute of Information Science, Academia Sinica, Taiwan.

# Introduction
RNA-Seq technology can provide high resolution measurement of expression and high sensitivity in detecting low abundance transcripts. However, RNA-seq data requires a huge amount of computational efforts since this technology can produce sequence reads on the order of million/billion base-pairs in a single day. The very first step is to align each sequence fragment against the reference genome.

DART adopts a divide-and-conquer strategy to handle RNA-Seq transcript alignments. Unlike most of read aligners that try to ex-tend a seed in both directions with a dynamic programming step, DART divides a read sequence into one or more segments to re-place the seed extension step. The experiment results on synthetic datasets and real datasets show that DART is a highly efficient aligner that yields the highest sensitivity and accuracy and spends the least amount of time among the selected aligners.

# Download

Please use the command 
  ```
  $ git clone https://github.com/hsinnan75/Dart.git
  ```
to download the package of DART.

# Installation

We provide the executable file, please type 

  ```
  $ ./dart index | aln [options]
  ```
to run the program. Or you can type 'make' to build the executable file.

# Usage

For indexing a reference genome, DART requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

  ```
  $ ./dart index -p ecoli Ecoli.fa
  ```

The above command is to index the genome file Ecoli.fa and store the index files begining with ecoli.

For mapping short reads, DART requires the the index files of the reference genome and at least one read file (two read files for the separated paired-end reads). Users should use -i to specify the prefix of the index files (including the directory path).

  ```
  $ ./dart aln -i ecoli -f ReadFile1.fq -f2 ReadFile2.fq -o out.sam
  ```

The above command is to run DART to align the paired-end reads in ReadFile1.fq and ReadFile2.fq with index files of ecoli. The output is redirected to out.sam.

# File formats

- a.Reference genome files

    All reference genome files should be in FASTA format.

- b.Read files

    All reads files should be in FASTA or FASTQ format. Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. The alignment result will not be different in either format.

    If paired-end reads are separated into two files, use -f and -f2 to indicate the two filenames. The i-th reads in the two files are paired. If paired-end reads are in the same file, use -p. The first and second reads are paired, the third and fourth reads are paired, and so on. For the latter case, use -p to indicate the input file contains paired-end reads.

- c.Output files

    Output is in standard SAM format. For reads aligned with reverse strand of reference genome, they are converted into obverse strand. More detailed information about SAM format, please refer to the SAMtools documents.
    
    We also output the predicted splice junctions (junctions.tab as the default filename or you may specify a filename with -j)

# Parameter setting

 ```

-t INT number of threads [16]

-i STR index prefix [BWT based (BWA), required]

-f STR read filename [required, fasta or fastq]

-f2 STR read filename2 [optional, fasta or fastq], f and f2 are files with paired reads

-p the input read file consists of interleaved paired-end sequences

-o STR alignment output

-j STR predicted splice junction filename [junctions.tab]

-m output multiple alignments

-unique output unique alignments

  ```
