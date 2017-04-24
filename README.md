DART: a fast and accurate RNA-seq mapper with divide and conquer strategy

Developers: Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu Institute of Information Science, Academia Sinica, Taiwan.

Introduction
DART is a read aligner for NGS RNA-Seq data developed by Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu. DART supports single-end and paired-end reads and multi-thread alignments. We describe the installation and Running instructions of DART below.

Instructions
1.Installation

We provide the executable file, please type './dart' to run the program. Or you can type 'make' to build the executable file.

2.Usage

For indexing a reference genome, DART requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

Ex. #./dart index -p ecoli Ecoli.fa

The above command is to index the genome file Ecoli.fa and store the index files begining with ecoli.

For mapping short reads, DART requires the the index files of the reference genome and at least one read file (two read files for the separated paired-end reads). Users should use -i to specify the prefix of the index files (including the directory path).

Ex. #./dart aln -i ecoli -f ReadFile1.fq -f2 ReadFile2.fq -o out.sam

The above command is to run DART to align the paired-end reads in ReadFile1.fq and ReadFile2.fq with index files of ecoli. The output is redirected to out.sam.

3.File formats

a.Reference genome files

All reference genome files should be in FASTA format.

b.Read files

All reads files should be in FASTA or FASTQ format. Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. The alignment result will not be different in either format.

If paired-end reads are separated into two files, use -f and -f2 to indicate the two filenames. The i-th reads in the two files are paired. If paired-end reads are in the same file, use -p. The first and second reads are paired, the third and fourth reads are paired, and so on. For the latter case, use -p to indicate the input file contains paired-end reads.

c.Output file

Output is in standard SAM format. For reads aligned with reverse strand of reference genome, they are converted into obverse strand. More detailed information about SAM format, please refer to the SAMtools documents.

Parameter setting

-t INT number of threads [16]

-i STR index prefix [BWT based (BWA), required]

-f STR read filename [required, fasta or fastq]

-f2 STR read filename2 [optional, fasta or fastq], f and f2 are files with paired reads

-p the input read file consists of interleaved paired-end sequences

-o STR alignment output

-m output multiple alignments

-unique output unique alignments
