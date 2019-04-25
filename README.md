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

# Changes
version 1.3.6: Discard the thread limit.

version 1.3.5: Add BAM format output.

version 1.3.4: Fix a bug on single-end mapping.

version 1.3.3: Fix a bug on paired-end mapping.

version 1.3.2: Replaced 0 exit codes with 1 and the corresponding 'Warning' with 'Error' for cases where program termination is not the expected result (revised by Rad Suchecki).

version 1.3.1: Add a "-v" option to show version of DART.

version 1.3.0: Fix a bug on reading input.

version 1.2.9: Fix a bug in SAM output.

version 1.2.8: Fix the alignment of segment pairs with poor sequence identity.

version 1.2.7: Add an update command.

version 1.2.6: Fix the alignment of segment pairs with multiple mismatches.

version 1.2.5: Fix the alignment when DNA sequences are shown in lower case.

version 1.2.4: Allow multiple read files as the input.

version 1.2.3: fix the bug when read number exceeds 2^23.

version 1.2.0: Add ksw2 and edlib alignment method to replace the Needleman-Wunsch algorithm.

version 1.1.2: fix a bug in the alignment report.

# Get updates

We update DART from time to time, please check if new version is available by using the following commands.

with Kart version 1.2.7 up
  ```
  $ ./dart update 
  ```
or
  ```
  $ git fetch
  $ git merge origin/master master
  $ make
  ```
# Compiling

To compile dart and the index tool, please change to dart's folder and just type 'make' to compile dart and bwt_index. If the compilation or the program fails, please contact me (arith@iis.sinica.edu.tw).

# Installation

We provide the executable file, please type 

  ```
  $ ./dart [options]
  ```
to run the program. Or you can type 'make' to build the executable file.

# Usage

To index a reference genome, DART requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

  ```
  $ ./bwt_index ref_file[ex.ecoli.fa] index_prefix[ex. Ecoli]
  ```
The above command is to index the genome file Ecoli.fa and store the index files begining with ecoli.

Please note that if you find bwt_index does not work in your computer system, you may also use bwa (http://bio-bwa.sourceforge.net/) to build the index files.
  ```
  $ ./bwa index -p index_prefix ref.fa
  ```

To map short reads, DART requires the the index files of the reference genome and at least one read file (two read files for the separated paired-end reads). Users should use -i to specify the prefix of the index files (including the directory path).

 case 1: standard sam output
  ```
 $ ./dart -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -o out.sam
  ```

 case 2: multiple input 
  ```
 $ ./dart -i ecoli -f ReadFileA_1.fq ReadFileB_1.fq ReadFileC_1.fq -f2 ReadFileA_2.fq ReadFileB_2.fq ReadFileC_2.fq -o out.sam
  ```

 case 3: bam output
  ```
 $ ./dart -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -bo out.bam
  ```

The above commands are to run DART to align the paired-end reads in ReadFile1.fq and ReadFile2.fq with index files of ecoli.

# File formats

- Reference genome files

    All reference genome files should be in FASTA format.

- Read files

    All reads files should be in FASTA or FASTQ format. FASTQ files can be compressed with gzip format. We do not support FASTA files with gzip compression.
    Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. The alignment result will not be different in either format.

    If paired-end reads are separated into two files, use -f and -f2 to indicate the two filenames. The i-th reads in the two files are paired. If paired-end reads are in the same file, use -p. The first and second reads are paired, the third and fourth reads are paired, and so on. For the latter case, use -p to indicate the input file contains paired-end reads.

- Output files

    Output is in standard SAM/BAM format. For reads aligned with reverse strand of reference genome, they are converted into obverse strand. More detailed information about SAM/BAM format, please refer to the SAMtools documents.

    We also output the predicted splice junctions (default: junctions.tab, or you may specify a filename with -j argument).

# Parameter setting

 ```

-t INT number of threads [16]

-i STR index prefix [BWT based (BWA), required]

-f STR read filename [required, fasta or fastq]

-f2 STR read filename2 [optional, fasta or fastq], f and f2 are files with paired reads

-p the input read file consists of interleaved paired-end sequences

-o STR alignment output [SAM]

-bo STR alignment output [BAM]

-j STR predicted splice junction filename [junctions.tab]

-m output multiple alignments

-intron INT the maximal intron size [500000]

-unique output unique alignments

  ```
