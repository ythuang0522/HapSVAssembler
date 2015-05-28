# Required thirdparty package
The bamtools are required to manually download while other programs are included.

# Compile
Type make in the src direcotry

# Execution
In the package folder, type ./Run_HapSV, and you will see the list of parameters

Version 1.1: released on May 22th, 2015

Usage: ./RunHapSV.sh -b <DIR> [OPTION]
  -f <int>  fasta of assembled contigs
  -1 <int>  fastq of 1st reads
  -2 <int>  fastq of 2nd reads
  -z <int>  insert size of paired-end reads
  -v <int>  standard deviation of insert size
  -l <int>  read length
  -h        help

For instance, given a contig file named "contig.fa" and fastq files of paried-end reads "R1.fq" and "R2.fq", simply type
./Run_HapSV -f contig.fa -1 R1.fq -2 R2.fq
