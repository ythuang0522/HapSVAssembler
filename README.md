# Required thirdparty package
The [bamtools][1] are required to manually download while other programs are included.

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

# Output
The output files are stored in a folder named haplotype_result. If there are multiple chromosomes/contigs, each of them is stored in a separate subfolder. Within each chromosome/contig subfolder, there are multiple output files explained below.

## haplotype.log
The main output file similar to the standard vcf format yet including block boundaries and SV information. Each row stands for a SNP or SV. The haplotype blocks are represented with a dash line. The following example shows three haplotype blocks composed of 6 SNPs and one SV insertion. Each row contains the chr./contig ID, variation type (SNP, Ins, Del, Inv), position, size (of SV if it is), allele on 1st haplotype, allele on 2nd haplotype.

    #ref_name	type	pos	size	hap0	hap1
    contig1001	SNP	396	-	G	G
    contig1001	SNP	405	-	T	T
    contig1001	SNP	470	-	G	G
    -
    contig1001	SNP	5398	-	C	T
    contig1001	Ins	5494	109	-	V
    contig1001	SNP	5565	-	G	A
    -
    contig1001	SNP	6380	-	G	A

## haplotype_variation.log
This output file stores the simplifed representation of paternal/maternal haplotypes sequences over SNPs only. There are two lines storing the paternal and maternal haplotype sequences, where the haplotype block boundary is represented using vertical line. e.g., three haplotype blocks composed of 22 SNPs.

    GTG|TCTCGTAATCACGATAGT|G
    GTG|AGATAACGCGGTGATAGT|A

## snp_site.log
The loci of detected SNPs.

## detected_sv.log
The loci of detected SVs.

## snp_matrix.log
The input SNP/SV matrix used in the constrained MEC model.

[1]: https://github.com/pezmaster31/bamtools
