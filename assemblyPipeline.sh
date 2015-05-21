#/bin/bash

# pair-end assembly (Kmers) -> amos merge -> scaffolding -> gapclosure
#=========================================
# parameter_setting
workDir=""
fq1=""
fq2=""
readLength=80
startKmer=31
endKmer=41
insertSize=250
variance=25
thread=1
# ========================================
# parse_arguments
while getopts 1:2:b:l:z:v:t:e:s:dirh OPTION
do
   case $OPTION in
   1) fq1="$OPTARG" ;;
   2) fq2="$OPTARG" ;;
   b) workDir="$OPTARG" ;;
   s) startKmer="$OPTARG" ;;
   e) endKmer="$OPTARG" ;;
   l) readLength="$OPTARG" ;;
   t) thread="$OPTARG" ;;
   z) insertSize="$OPTARG" ;;
   v) variance="$OPTARG" ;;
   h) echo ""
      echo "Version 1.0: released on November 13th, 2011"
      echo ""
      echo "Usage: $0 -b <DIR> [OPTION]"
      echo "  -b <DIR>  the path of HapSVAssembler package"
      echo "  -1 <PATH/filename> first end of paired-end reads in FASTQ format"
      echo "  -2 <PATH/filename> second end of paired-end reads in FASTQ format"
      echo "  -l <int>  number of known nucleotides on each ends"
      echo "  -z <int>  insert size of paired-end reads"
      echo "  -v <int>  standard deviation of insert size"
      echo "  -s <int>  start kmer size should be odd"
      echo "  -e <int>  end kmer size should be odd"
      echo "  -t <int>  number of threads"
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
      echo "Usage: $0 -b <DIR> [OPTION]"
      echo "  -b <DIR>  the path of HapSVAssembler package"
      echo "  -1 <PATH/filename> first end of paired-end reads in FASTQ format"
      echo "  -2 <PATH/filename> second end of paired-end reads in FASTQ format"
      echo "  -l <int>  number of known nucleotides on each ends"
      echo "  -z <int>  insert size of paired-end reads"
      echo "  -v <int>  standard deviation of insert size"
      echo "  -s <int>  start kmer size"
      echo "  -e <int>  end kmer size"
      echo "  -t <int>  number of threads"
      echo "  -h        help"
      echo ""
      exit 0
fi

#==========================================
# software path
SOAPtool=$workDir/bin/SOAPdenovo127mer
toAmos=$workDir/bin/toAmos
minimus2=$workDir/bin/minimus2
combineTwoSOAP=$workDir/bin/combineTwoSOAP
Filter=$workDir/bin/filterContigLength
BWA=$workDir/bin/bwa
N50=$workDir/bin/n50.sh
Prepare=$workDir/bin/prepare
Gapcloser=$workDir/bin/GapCloser
#==========================================


cd $workDir/denovo_result

#ln -s $fq1 $workDir/denovo_result/testGenome_1.fastq
#ln -s $fq2 $workDir/denovo_result/testGenome_2.fastq
#ln -s $fq1 $workDir/genome_result/testGenome_1.fastq
#ln -s $fq2 $workDir/genome_result/testGenome_2.fastq

######################################
# create contig config file
        echo -n "" > config_contig
        echo "max_rd_len=$readLength" >> config_contig
        echo "" >> config_contig
        echo "[LIB]"  >> config_contig
        echo "asm_flags=1"  >> config_contig
        echo "q1=$workDir/$fq1" >> config_contig
        echo "q2=$workDir/$fq2" >> config_contig

######################################
# create scaffoding config file
	echo -n "" > config_scaf
	echo "max_rd_len=$readLength" > config_scaf
	echo "" >> config_scaf
	echo "[LIB]"  >> config_scaf
	echo "avg_ins=$insertSize" >> config_scaf
	echo "asm_flags=3"  >> config_scaf
	echo "q1=$workDir/$fq1" >> config_scaf
	echo "q2=$workDir/$fq2" >> config_scaf

######################################
# create gap close config file
	echo -n "" > config_close
	echo "max_rd_len=$readLength" > config_close
	echo "" >> config_close
	echo "[LIB]"  >> config_close
	echo "avg_ins=$insertSize" >> config_close
	echo "asm_flags=4"  >> config_close
	echo "q1=$workDir/$fq1" >> config_close
	echo "q2=$workDir/$fq2" >> config_close


# assembly
echo -n "" > ctemp.fa
for ((i=$startKmer; i<=$endKmer; i=i+2 ))
do
        mkdir Kmer$i
        cd Kmer$i

# start denovo assembly by soap tool
        $SOAPtool all -K $i -p $thread -s ../config_contig -o Kmer$i

# remove non-important files
        rm *.Arc *.ContigIndex *.edge *.gapSeq *.links *.newContigIndex *.peGrads *.preArc *.preGraphBasic *.readInGap *.readOnContig *.scaf *.scaf_gap *.vertex

# keep greater than 500bp
        $Filter Kmer$i.contig out 500
        mv out_backbone_in.fasta g500_Kmer$i.contig
        rm out_in.sanger.fasta

# create ctemp.fa
	cat g500_Kmer$i.contig >> ../ctemp.fa

	cd ..

#$N50 Kmer$i/g500_Kmer$i.contig
done

#$N50 ctemp.fa
######################################
# merge different kmer contigs by amos
$combineTwoSOAP ctemp.fa
$toAmos -s ctemp.fa.reorder -o ctemp.afg
$minimus2 -D OVERLAP=80 ctemp
#rm *.contig *.coords *.delta *.fa *.ovl *.OVL *.ref.seq *.qry.seq *.log *.singletons
#$N50 ctemp.fasta

######################################
# scaffolding
$Prepare -g peScaf -K $startKmer -c ctemp.fasta
$SOAPtool map -s config_scaf -p $thread -g peScaf
$SOAPtool scaff -g peScaf -p $thread
#$N50 peScaf.scafSeq

######################################
# gap closure
$Gapcloser -b config_close -o peCloser -a peScaf.scafSeq
#$N50 peCloser

rm -r ctemp.fa ctemp.contig *.coords *.delta *.ovl *.OVL *.qry.seq *.ref.seq *.runAmos.log *.bnk *.singletons *.afg
rm *.fill *.Arc *.ContigIndex *.edge *.gapSeq *.links *.newContigIndex *.peGrads *.preArc *.preGraphBasic *.readInGap *.readOnContig *.scaf *.scaf_gap *.updated.edge *.vertex *.newContigIndex *.path *.markOnEdge *.conver

