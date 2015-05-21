workDir=`pwd`
fasta=""
fq1=""
fq2=""

readLength=100
insertSize=200
variance=50

while getopts b:f:1:2:z:l:v:dirh OPTION
do
   case $OPTION in
    b) workDir="$OPTARG" ;;
    f) fasta="$OPTARG" ;;
    1) fq1="$OPTARG" ;;
    2) fq2="$OPTARG" ;;
    l) readLength="$OPTARG" ;;
    z) insertSize="$OPTARG" ;;
    v) variance="$OPTARG" ;;
    h) echo ""
       echo "Version 1.1: released on 5/22, 2015"
       echo "Written by YTH and SYC"
       echo "Usage: $0 -b <DIR> [OPTION]"
       echo "  -f <int>  fasta of assembled contigs"  
       echo "  -1 <int>  fastq of 1st reads"  
       echo "  -2 <int>  fastq of 2nd reads"  
       echo "  -z <int>  insert size of paired-end reads"  
       echo "  -v <int>  standard deviation of insert size"  
       echo "  -l <int>  read length"  
       echo "  -h        help"
       echo ""
       exit 0
       ;;
   esac
done

if [ $OPTIND -lt 4 ]; then
      echo ""
      echo "Version 1.1: released on May 22th, 2015"
      echo ""
      echo "Usage: $0 -b <DIR> [OPTION]"
       echo "  -f <int>  fasta of assembled contigs"  
       echo "  -1 <int>  fastq of 1st reads"  
       echo "  -2 <int>  fastq of 2nd reads"  
      echo "  -z <int>  insert size of paired-end reads"
      echo "  -v <int>  standard deviation of insert size"
      echo "  -l <int>  number of known nucleotides on each ends"
        
      echo "  -h        help"
      echo ""
      exit 0
fi

#=====bwa========================================================================================================================
before=`date '+%s'`

rm -f $workDir/bwa_result/reference.fasta $workDir/bwa_result/testGenome_1.fastq $workDir/bwa_result/testGenome_2.fastq
#ln -s $fasta $workDir/bwa_result/reference.fasta
#ln -s $fq1 $workDir/bwa_result/testGenome_1.fastq
#ln -s $fq2 $workDir/bwa_result/testGenome_2.fastq

#cd $workDir/bwa_result/
rm -f $workDir/bwa_result/*.bam $workDir/bwa_result/*.sam

#fold -w 100 reference.fasta > reference_fold.fasta
#$workDir/bin/bwa index reference_fold.fasta
#$workDir/bin/bwa mem reference_fold.fasta testGenome_1.fastq testGenome_2.fastq -t 31 > testGenome.sam
$workDir/bin/bwa index $fasta 
$workDir/bin/bwa mem $fasta $fq1 $fq2 -t 31 > $workDir/bwa_result/testGenome.sam

after=`date '+%s'`
elapsed=`echo "$after - $before" | bc`
hours=`echo "elapsed / 3600" | bc`
elapsed=`echo "$elapsed - $hours * 3600" | bc`
minutes=`echo "$elapsed / 60" | bc`
seconds=`echo "$elapsed - $minutes * 60" | bc`
if [ $hours != 0 ]; then
	echo "[bwa] overall time elapses: $hours Hours $minutes Minutes $seconds Seconds";
elif [ $minutes != 0 ]; then
	echo "[bwa] overall time elapses: $minutes Minutes $seconds Seconds";
else
	echo "[bwa] overall time elapses: $seconds Seconds";
fi
echo ""
#================================================================================================================================


#=====samtools===================================================================================================================
cd $workDir/bwa_result/
before=`date '+%s'`
#echo "[samtools_faidx] indexing fasta..."
#samtools faidx reference_fold.fasta
echo "[samtools_view] converting sam to bam..."
printf " -"
samtools view -Sb testGenome.sam > testGenome.bam
echo "[samtools_sort] sorting bam..."
samtools sort testGenome.bam testGenome_sorted

after=`date '+%s'`
elapsed=`echo "$after - $before" | bc`
hours=`echo "elapsed / 3600" | bc`
elapsed=`echo "$elapsed - $hours * 3600" | bc`
minutes=`echo "$elapsed / 60" | bc`
seconds=`echo "$elapsed - $minutes * 60" | bc`
if [ $hours != 0 ]; then
        echo "[samtools] overall time elapses: $hours Hours $minutes Minutes $seconds Seconds";
	elif [ $minutes != 0 ]; then
	        echo "[samtools] overall time elapses: $minutes Minutes $seconds Seconds";
	else
	        echo "[samtools] overall time elapses: $seconds Seconds";
fi
echo ""
#================================================================================================================================

#================================================================================================================================
cd $workDir/bwa_result/
$workDir/bin/bamtools/bin/bamtools split -in testGenome_sorted.bam -reference
$workDir/bin/bamtools/bin/bamtools split -in testGenome.bam -reference
#================================================================================================================================



cd $workDir
ContigIDs=($(grep ">" "$fasta" | tr -d '>' | cut -f 1 -d ' '))
for contigID in "${ContigIDs[@]}"
do
	echo "Processing " $contigID "\n"
	grep "$contigID" "$fasta" -A 1 > tmp_contig.fasta
    $workDir/bin/HapSVAssembler2.sh -b $workDir -1 testGenome.REF_$contigID.bam -2 testGenome_sorted.REF_$contigID.bam -a $workDir/tmp_contig.fasta -z $insertSize -v $variance -l $readLength 
    rm $workDir/haplotype_result/paired-align.dist $workDir/haplotype_result/sam_summary.log tmp_contig.fasta*
    mkdir $workDir/haplotype_result/$contigID
    mv $workDir/haplotype_result/*.log $workDir/haplotype_result/$contigID
    rm $workDir/haplotype_result/snp.vcf
done
