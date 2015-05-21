#=====parameter_setting==========================================================================================================
workDir=""
fasta=""
fq1=""
fq2=""
genome_size=10000
read_length=80
start=0
insert_size=250
error_rate=0
snp_rate=100
variance=25
coverage=20
difference=1
fasta_name=""
del_flag=""
ins_flag=""
rev_flag=""
pop_size=10
generation=5
mutation_rate=80
TMPDIR="tmpDir"
#================================================================================================================================

#=====parse_arguments============================================================================================================
while getopts a:1:2:b:f:l:S:z:v:G:c:p:g:m:e:s:D:dirh OPTION
do
   case $OPTION in
   a) fasta="$OPTARG" ;;
   1) fq1="$OPTARG" ;;
   2) fq2="$OPTARG" ;;
   b) workDir="$OPTARG" ;;
   f) fasta_name="-f $OPTARG" ;;
   S) start="$OPTARG" ;;
   G) genome_size="$OPTARG" ;;
   l) read_length="$OPTARG" ;;
   z) insert_size="$OPTARG" ;;
   v) variance="$OPTARG" ;;
   c) coverage="$OPTARG" ;;
   s) snp_rate="$OPTARG" ;;
   D) difference="$OPTARG" ;;
   d) del_flag="-d" ;;
   i) ins_flag="-i" ;;
   r) rev_flag="-r" ;;
   e) error_rate="$OPTARG" ;;
   p) pop_size="$OPTARG" ;;
   g) generation="$OPTARG" ;;
   m) mutation_rate="$OPTARG" ;;
   h) echo ""
      echo "Version 1.0: released on November 13th, 2011"
      echo ""
      echo "Usage: $0 -b <DIR> [OPTION]" 
      echo "  -b <DIR>  the path of HapSVAssembler package"
      echo "  -a <FILE> reference genome/chromosome in FASTA format"
      echo "  -1 <FILE> first end of paired-end reads in FASTQ format"
      echo "  -2 <FILE> second end of paired-end reads in FASTQ format"
      echo ""
      echo "  [OPTION]"
      echo "  -f <FILE> filename of reference genome in FASTA format"
      echo "  -S <int>  started-cutting position on reference [Default 0] (required only when -f <FILE> is set)"
      echo "  -G <int>  genome size for simulated data [Default 10000]"
      echo "  -l <int>  number of known nucleotides on each ends [Default 80]"
      echo "  -z <int>  insert size of paired-end reads [Default 250]"
      echo "  -v <int>  standard deviation of insert size [Default 25]"
      echo "  -c <int>  coverage [Default 20]"
      echo "  -s <int>  SNP rate (1/s, s=[100,10000]) in diploid genome [Default 100]"
      echo "  -D <int>  difference rate (D%, s=[1,5]) between two homologous chromosomes [Default 1]"
      echo "  -d        create deletion polymorphism [Default NO]"
      echo "  -i        create insertion polymorphism [Default NO]"
      echo "  -r        create reversion polymorphism [Default NO]"
      echo ""
      echo "  -e <int>  error rate (e%, e=[0,100]) in SNP matrix [Default 0]"
      echo ""
      echo "  -p <int>  population size for GA-pipeline [Default 10]"
      echo "  -g <int>  generation for GA-pipeline [Default 5]"
      echo "  -m <int>  mutation rate (m%, m=[0,100]) for GA-pipeline [Default 80]"
      echo "  -h        help"
      echo ""
      exit 0
      ;;
   esac
done

if [ $OPTIND -lt 3 ]; then 
	echo ""
        echo "Version 1.0: released on November 13th, 2011"
	echo ""
	echo "Usage: $0 -b<DIR> [OPTION]" 
	echo "  -b <DIR> the path of HapSVAssembler package"
	echo ""
	echo "  [OPTION]"
	echo "  -f <file> filename of reference genome in FASTA format"
	echo "  -G <int>  genome size for simulated data [Default 10000]"
	echo "  -S <int>  started-cutting position on reference [Default 0] (required only when -f <FILE> is set)"
	echo "  -l <int>  number of known nucleotides on each ends [Default 80]"
	echo "  -z <int>  insert size of paired-end reads [Default 250]"
	echo "  -v <int>  standard deviation of insert size [Default 25]"
	echo "  -c <int>  coverage [Default 20]"
	echo "  -s <int>  SNP rate (1/s, s=[100,10000]) in diploid genome [Default 100]"
	echo "  -D <int>  difference rate (D%, s=[1,5]) between two homologous chromosomes [Default 1]"
	echo "  -d        create deletion polymorphism [Default NO]"
	echo "  -i        create insertion polymorphism [Default NO]"
	echo "  -r        create reversion polymorphism [Default NO]"
	echo ""
	echo "  -e <int>  e% error rate in SNP matrix [Default 0]"
	echo ""
	echo "  -p <int>  population size for GA-pipeline [Default 10]"
	echo "  -g <int>  generation for GA-pipeline [Default 5]"
	echo "  -m <int>  mutation rate(m/100) for GA-pipeline [Default 80]"
	echo "  -h        help"
	echo ""
	exit 0
fi
#================================================================================================================================

#=====genomeCreator==============================================================================================================
#echo ""
#cd $workDir/genome_result/
#$workDir/bin/genomeCreator $fasta_name -s $start -g $genome_size -c $coverage -l $read_length -s $snp_rate -D $difference $del_flag $ins_flag $rev_flag > tmp
#if [ $? != 0 ]; then
#	rm tmp
#	exit 0
#fi
#cat tmp
#read1=`grep "chromatid1" tmp | cut -d ':' -f 2`
#read2=`grep "chromatid2" tmp | cut -d ':' -f 2`
#rm tmp
#================================================================================================================================

#=====wgsim======================================================================================================================
#before=`date '+%s'`
#$workDir/bin/wgsim reference.fasta tmp1_1.fq tmp1_2.fq -e 0 -d $insert_size -s $variance -N $read1 -1 $read_length \
#      -2 $read_length -r 0 -R 0 -X 0
#$workDir/bin/wgsim reference_chromatid2.fasta tmp2_1.fq tmp2_2.fq -e 0 -d $insert_size -s $variance -N $read1 \
#      -1 $read_length -2 $read_length -r 0 -R 0 -X 0
#cat tmp1_1.fq tmp2_1.fq > testGenome_1.fastq
#cat tmp1_2.fq tmp2_2.fq > testGenome_2.fastq
#rm tmp*
#after=`date '+%s'`
#elapsed=`echo "$after - $before" | bc`
#hours=`echo "elapsed / 3600" | bc`
#elapsed=`echo "$elapsed - $hours * 3600" | bc`
#minutes=`echo "$elapsed / 60" | bc`
#seconds=`echo "$elapsed - $minutes * 60" | bc`
#if [ $hours != 0 ]; then
#	echo "[wgsim] time elapses: $hours Hours $minutes Minutes $seconds Seconds";
#elif [ $minutes != 0 ]; then
#	echo "[wgsim] time elapses: $minutes Minutes $seconds Seconds";
#else
#	echo "[wgsim] time elapses: $seconds Seconds";
#fi
#echo "[wgsim]$read1 paired-end reads from chromatid1 have been generated."
#echo "[wgsim]$read2 paired-end reads from chromatid2 have been generated."
#echo ""
#================================================================================================================================


#=====preprocess=================================================================================================================
export PATH=$PATH:$workDir/bin/
#rm -f reference.fasta testGenome_1.fastq testGenome_2.fastq
#ln -s $fasta reference.fasta
#ln -s $fq1 testGenome_1.fastq
#ln -s $fq2 testGenome_2.fastq
#================================================================================================================================


#=====bwa========================================================================================================================
before=`date '+%s'`
cd $workDir/bwa_result/
rm -f reference.fasta testGenome_1.fastq testGenome_2.fastq
ln -s $workDir/genome_result/reference.fasta
#ln -s $workDir/denovo_result/result.contig ./reference.fasta
ln -s $workDir/genome_result/testGenome_1.fastq
ln -s $workDir/genome_result/testGenome_2.fastq
fold -w 100 reference.fasta > reference_fold.fasta
$workDir/bin/bwa index -a is reference_fold.fasta
$workDir/bin/bwa aln reference_fold.fasta testGenome_1.fastq -f testGenome_1.sai -t 31
$workDir/bin/bwa aln reference_fold.fasta testGenome_2.fastq -f testGenome_2.sai -t 31
$workDir/bin/bwa sampe -a 1000 reference_fold.fasta testGenome_1.sai testGenome_2.sai testGenome_1.fastq testGenome_2.fastq -f testGenome.sam
$workDir/bin/bwa samse reference_fold.fasta testGenome_1.sai testGenome_1.fastq -f testGenome_1.sam
$workDir/bin/bwa samse reference_fold.fasta testGenome_2.sai testGenome_2.fastq -f testGenome_2.sam
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
before=`date '+%s'`
echo "[samtools_faidx] indexing fasta..."
samtools faidx reference_fold.fasta
echo "[samtools_view] converting sam to bam..."
printf " -"
samtools view -Sb testGenome_1.sam > testGenome_1.bam
printf " -"
samtools view -Sb testGenome_2.sam > testGenome_2.bam
echo "[samtools_sort] sorting bam..."
samtools sort testGenome_1.bam testGenome_1_sorted
samtools sort testGenome_2.bam testGenome_2_sorted
echo "[samtools_mpileup] calling SNPs/INDELs..."
printf " -"
samtools mpileup -uf reference_fold.fasta testGenome_1_sorted.bam testGenome_2_sorted.bam | bcftools view -vcg - > snp.vcf
#grep -v "#" snp.vcf | cut -f 1,2,4,5 | awk '{printf "%s\t%d\t%c\t%c\n",$1,$2-1,$3,$4}' > snp_site.log
grep -v "#" snp.vcf | cut -f 1,2,4,5 | awk '($3=="A"||$3=="C"||$3=="G"||$3=="T")&&($4=="A"||$4=="C"||$4=="G"||$4=="T")' | awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4}' > snp_site.log
#grep -v "#" snp.vcf | cut -f 1,2,4,5 | awk '{print $1 "\t" $2-1 "\t" $3 "\t" $4}' > snp_site.log
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

#=====svDetector=================================================================================================================
cd $workDir/haplotype_result/
echo "[preprocessing_sam] converting SAM..."
awk '$1!="@PG" && $1!="@SQ"'  $workDir/bwa_result/testGenome.sam | cut -f 1,2,4,6,8,9,10 > sam_summary.log
awk 'BEGIN{OFS="\t"} {print $4,$2,$6,$3,$5}' sam_summary.log > paired-align.dist
r=`wc -l paired-align.dist | cut -d ' ' -f 1`
cut -f 1,4,6 $workDir/bwa_result/testGenome.sam | sort -T $workDir/$TMPDIR -k 2,2 -n | awk '$2!= 0' | awk '$3!="*"' | grep 'S' | grep -v '@' > tmp
grep -G '[0-9]*S[0-9]' tmp | grep -v -G '[0-9]*S[0-9]*M[0-9]' | cut -f 2 | uniq -c | awk '$1>=2' | awk '{print $2-1"\t"1}' > rightBP
cut -f 2,3 tmp | sed 's/M/\t/g' | awk '$3!=""' | awk '$4==""' | grep -v -G '[0-9]*S[0-9]' | awk '{print $1+$2-1}' | sort -n | uniq -c | awk '$1>=2' | awk '{print $2"\t"2}' > leftBP
cat leftBP rightBP | sort -n > breakPoint.log
b=`wc -l breakPoint.log | cut -d ' ' -f 1`
rm tmp *BP
$workDir/bin/svDetector -z $insert_size -v $variance -l $read_length -r $r -b $b
echo ""
#================================================================================================================================

#=====snpMatrixCreator===========================================================================================================
rm -f snp_site.log
cp $workDir/bwa_result/snp_site.log .
cp $workDir/bwa_result/snp.vcf .
sort -T $workDir/$TMPDIR -n -k 3,3 Detected_SV_With_Count.log | cut -f 1-4 > detected_sv.log
s=`wc -l snp_site.log | cut -d ' ' -f 1`
r=`wc -l sam_summary.log | cut -d ' ' -f 1`
V=`wc -l detected_sv.log | cut -d ' ' -f 1`
$workDir/bin/snpMatrixCreator -s $s -r $r -V $V -z $insert_size -v $variance -l $read_length
sort -T $workDir/$TMPDIR -n -k 4,4 tmpMatrix > snp_matrix.log
rm tmpMatrix
echo ""
#================================================================================================================================

#=====hapAssembler===============================================================================================================
r=`wc -l snp_matrix.log | cut -d ' ' -f 1`
c=`wc -l cluster.log | cut -d ' ' -f 1`
$workDir/bin/hapAssembler -s $s -r $r -c $c -V $V -p $pop_size -g $generation -m $mutation_rate
#================================================================================================================================
