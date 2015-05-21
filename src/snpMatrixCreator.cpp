#include <unistd.h>
#include <fstream>
#include <string.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "SNP.h"
#include "SV.h"

using namespace std;

int snpCount = 0; 
int readCount = 0;
int svCount = 0;
int insertSize = 250;
int readLength = 75;
int variance = 25;
int errorRate = 0;
int clusterIndex = 0; //overall clusters in SNP matrix
int validIndex = 0;   //number of reads that at least one-end matches with reference genome perfectly
int rightAlignNum = 0; //number of reads that only left-end matches with reference genome perfectly
int offset = 10; //offset of indel
struct READ
{
   int alignType;  //1: only second end matching, 2: only first end matching, 3: matching both ends
   int pos[2];     //[first end posistion,second end position]
   int mappingDist;
   char *seq[2];   //[first end sequences,second end sequences]
};
SNP *snpSite; //dynamic SNP array
SV *svSite;   //dynamic SV array
int *snpSet;  //the cluster number of corresponding SNP site(e.g., snpSet[10]=3 means SNP10 belongs to CLUSTER3)
struct READ *readList;
int **snpMatrix;

void readSNP(); 
void readSV();
void readSAM();
void sortRead();
void initialReadStruct();
void createSNPmatrix();
bool between(int start, int leftMost, int rightMost, int *startSNP, int *endSNP);
static int usage();
int comparePos1(const void *s1, const void *s2) { return ((*(struct READ*)s1).pos[0] - (*(struct READ*)s2).pos[0]); }
int comparePos2(const void *s1, const void *s2) { return ((*(struct READ*)s1).pos[1] - (*(struct READ*)s2).pos[1]); }
void printMatrix();
void printCluster();

int main(int argc, char *argv[])
{
   double startT, endT, overallTime = 0, sectionTime = 0;
   char c;

   while((c = getopt(argc,argv,"s:r:V:z:v:l:h")) >= 0)
   {
      switch(c)
      {
         case 's': snpCount = atoi(optarg); break;
         case 'r': readCount = atoi(optarg); break;
	 case 'V': svCount = atoi(optarg); break;
	 case 'z': insertSize = atoi(optarg); break;
	 case 'v': variance = atoi(optarg); break;
	 case 'l': readLength = atoi(optarg); break;
         case 'h': break;
         default : cout << "Unrecognized option: -" << c << endl; return 1;
      }
   }
   if(optind < 13) return usage();

   snpSite = new SNP [snpCount];
   svSite = new SV [svCount];
   readList = new struct READ [readCount/2];
   snpSet = new int [snpCount];
   for(int i=0;i<snpCount;i++) snpSet[i] = -1;

   printf( "[snpMatrixCreator_arguments] %s -s %d -r %d -V %d ",argv[0],snpCount,readCount,svCount);
   printf( "-z %d -v %d -l %d \n",insertSize,variance,readLength);

   startT = clock();
   readSNP();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[snpMatrixCreator_readSNP] time elapses: %.2lf Seconds\n",sectionTime);
   cout << "[snpMatrixCreator_readSNP] " << readCount << " reads have been loaded." << endl;

   startT = clock();
   readSV();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[snpMatrixCreator_readSV] time elapses: %.2lf Seconds\n",sectionTime);
   cout << "[snpMatrixCreator_readSV] " << svCount << " detected structural variations have been loaded." << endl;

   startT = clock();
   readSAM();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[snpMatrixCreator_readSAM] time elapses: %.2lf Seconds\n",sectionTime);
   cout << "[snpMatrixCreator_readSAM] SAM file has been processed." << endl;

   startT = clock();
   qsort(readList,validIndex,sizeof(struct READ),comparePos1);
   qsort(readList,rightAlignNum,sizeof(struct READ),comparePos2);
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[snpMatrixCreator_sorting] time elapses: %.2lf Seconds\n",sectionTime);
   cout << "[snpMatrixCreator_sorting] read positions have been sorted." << endl;

   startT = clock();
   cout << "[snpMatrixCreator_createSNPmatrix] creating SNP matrix... ";
   cout.flush();
   createSNPmatrix();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   if(sectionTime/3600 > 1) 
   {
      printf( "%.0lf Hours ",sectionTime/3600);
      sectionTime = sectionTime - floor(sectionTime/3600)*3600;
   }
   if(sectionTime/60 > 1) 
   {
      printf( "%.0lf Minutes ",sectionTime/60);
      sectionTime = sectionTime - floor(sectionTime/60)*60;
   }
   printf( "%.2lf Seconds\n",sectionTime);

   startT = clock();
   printMatrix();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[snpMatrixCreator_printMatrix] time elapses: %.2lf Seconds\n",sectionTime);
   cout << "[snpMatrixCreator_printMatrix] SNP matrix has written into disk." << endl;

   startT = clock();
   printCluster();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[snpMatrixCreator_printCluster] time elapses: %.2lf Seconds\n",sectionTime);
   cout << "[snpMatrixCreator_printCluster] SNP cluster has written into disk." << endl;

   printf( "[snpMatrixCreator] overall time elapses: ");
   if(overallTime/3600 > 1)
   {
      printf( "%.0lf Hours ",overallTime/3600);
      overallTime = overallTime - floor(overallTime/3600)*3600;
   }
   if(overallTime/60 > 1)
   {
      printf( "%.0lf Minutes ",overallTime/60);
      overallTime = overallTime - floor(overallTime/60)*60;
   }
   printf( "%.2lf Seconds\n",overallTime);

   return 0;
}

static int usage()
{
   cout << "\nVersion 1.0: released on November 3rd, 2011" << endl;
   cout << "\nUsage: snpMatrixCreator -s <int> -r <int> -m <int> -z <int> -V <int> -l <int>" << endl;
   cout << "  -s <int> number of SNP sites" << endl;
   cout << "  -r <int> number of reads in SAM" << endl;
   cout << "  -V <int> number of structural variations" << endl;
   cout << "  -z <int> insert size of paired-end reads" << endl;
   cout << "  -v <int> standard deviation of insert size" << endl;
   cout << "  -l <int> number of known nucleotides on each ends" << endl;
   cout << "  -h       help" << endl << endl;
   return 1;
}

void readSNP()
{
   ifstream fin;
   char buf[1024], allele[2];
   int pos;

   fin.open("snp_site.log");
   for(int i=0;i<snpCount;i++)
   {
      fin.getline(buf,1024);
      strtok(buf,"\t");
      pos = atoi(strtok(NULL,"\t"));
      allele[0] = *(strtok(NULL,"\t"));
      allele[1] = *(strtok(NULL,"\n"));
      snpSite[i].insertSNP(pos,allele[0],allele[1]);
   }
   fin.close();
}

void readSV()
{
   ifstream fin;
   char buf[1024], *tmp;
   int type, size, pos[2];

   fin.open("detected_sv.log");
   for(int i=0;i<svCount;i++)
   {
      fin.getline(buf,1024);
      tmp = strtok(buf,"\t");
      if(strcmp(tmp,"Deletion") == 0) type = 1;
      else if(strcmp(tmp,"Insertion") == 0) type = 2;
      else if(strcmp(tmp,"Inversion") == 0) type = 3;
      size = atoi(strtok(NULL,"\t"));
      pos[0] = atoi(strtok(NULL,"\t"));
      pos[1] = atoi(strtok(NULL,"\n"));
      svSite[i].insertSV(pos[0],pos[1],size,type);
   }
   fin.close();
}

void readSAM()
{
   ifstream fin;
   char buf[1024], readName[2][64], cigar[2][16], seq[2][readLength+1], *tmpPtr, *testCigar, digit[3];
   int flag[2], selfPos[2], matePos[2], mappingDist[2];
   int i, j, index = 0;
   int tmpLength = readLength;

   //transforming numerals into separated-decile characters
   //tmpLength: original read length(e.g., 80)
   //digit: separated-decile array
   //ascii 0 is 48
   while(tmpLength/10 != 0)
   {
      digit[index] = tmpLength%10 + 48;
      tmpLength = tmpLength/10;
      index++;
   }
   digit[index] = tmpLength + 48;

   //generating perfect matching CIGAR
   //testCigar: perfect matching CIGAR(e.g, 80M)
   testCigar = new char [index+2];
   for(i=index;i>=0;i--) testCigar[index-i] = digit[i];
   testCigar[index+1] = 'M';
   testCigar[index+2] = '\0';

   fin.open("sam_summary.log");
   for(i=0;i<readCount;)
   {
      for(j=0;j<2;j++) //read pair-end information
      {
	 fin.getline(buf,1024);
	 i++;
	 strcpy(readName[j],strtok(buf,"\t"));
	 readName[j][strlen(readName[j])-1] = '\0';
	 flag[j] = atoi(strtok(NULL,"\t"));
	 selfPos[j] = atoi(strtok(NULL,"\t")) - 1;
	 strcpy(cigar[j],strtok(NULL,"\t"));
	 matePos[j] = atoi(strtok(NULL,"\t")) - 1;
	 mappingDist[j] = atoi(strtok(NULL,"\t"));
	 strcpy(seq[j],strtok(NULL,"\n"));
      }
      if(strcmp(readName[0],readName[1]) != 0)
      {
	 fin.getline(buf,1024);
	 i++;
	 continue;
      }
      if(strncmp(cigar[0],testCigar,index+1) != 0 && strncmp(cigar[1],testCigar,index+1) != 0) continue;
      if(abs(mappingDist[0]) > 5000) continue;
      if(strncmp(cigar[0],testCigar,index+1) == 0 && strncmp(cigar[1],testCigar,index+1) == 0) 
      {
	 initialReadStruct();
	 if((flag[0] == 99 && flag[1] == 147) || (flag[0] == 97 && flag[1] == 145)) //(+.-)
	 {
	    readList[validIndex].seq[0] = new char [readLength+1];
	    readList[validIndex].seq[1] = new char [readLength+1];
	    readList[validIndex].alignType = 3;
	    readList[validIndex].pos[0] = selfPos[0];
	    readList[validIndex].pos[1] = matePos[0];
	    readList[validIndex].mappingDist = mappingDist[0];
	    strcpy(readList[validIndex].seq[0],seq[0]);
	    strcpy(readList[validIndex].seq[1],seq[1]);
	    validIndex++;
	 }
	 else if((flag[0] == 83 && flag[1] == 163) || (flag[0] == 81 && flag[1] == 161)) //(-,+)
	 {
	    readList[validIndex].seq[0] = new char [readLength+1];
	    readList[validIndex].seq[1] = new char [readLength+1];
	    readList[validIndex].alignType = 3;
	    readList[validIndex].pos[0] = selfPos[1];
	    readList[validIndex].pos[1] = matePos[1];
	    readList[validIndex].mappingDist = mappingDist[1];
	    strcpy(readList[validIndex].seq[0],seq[1]);
	    strcpy(readList[validIndex].seq[1],seq[0]);
	    validIndex++;
	 }
	 else if(flag[0] == 65 && flag[1] == 129) //(+.+)
	 {
	    readList[validIndex].alignType = 2;
	    readList[validIndex].seq[0] = new char [readLength+1];
	    readList[validIndex].seq[1] = NULL;
	    if(mappingDist[0] > 0) //first end is on left-handed side
	    {
	       readList[validIndex].pos[0] = selfPos[0];
	       strcpy(readList[validIndex].seq[0],seq[0]);
	    }
	    else //first end is on right-handed side
	    {
	       readList[validIndex].pos[0] = selfPos[1];
	       strcpy(readList[validIndex].seq[0],seq[1]);
	    }
	    validIndex++;
	 }
	 else if(flag[0] == 113 && flag[1] == 177) //(-,-)
	 {
	    readList[validIndex].alignType = 1;
	    readList[validIndex].seq[1] = new char [readLength+1];
	    readList[validIndex].seq[0] = NULL;
	    if(mappingDist[0] < 0) //first end is on right-handed side
	    {
	       readList[validIndex].pos[1] = selfPos[0];
	       strcpy(readList[validIndex].seq[1],seq[0]);
	    }
	    else //first end is on left-handed side
	    {
	       readList[validIndex].pos[1] = selfPos[1];
	       strcpy(readList[validIndex].seq[1],seq[1]);
	    }
	    rightAlignNum++;
	    validIndex++;
	 }
      }
      else if(strncmp(cigar[0],testCigar,index+1) == 0)
      {
	 initialReadStruct();
	 if(flag[0] == 73 || flag[0] == 99)
	 {
	    readList[validIndex].seq[0] = new char [readLength+1];
	    readList[validIndex].seq[1] = NULL;
	    readList[validIndex].alignType = 2;
	    readList[validIndex].pos[0] = selfPos[0];
	    strcpy(readList[validIndex].seq[0],seq[0]);
	    validIndex++;
	 }
	 else if(flag[0] == 89 || flag[0] == 83)
	 {
	    readList[validIndex].seq[1] = new char [readLength+1];
	    readList[validIndex].seq[0] = NULL;
	    readList[validIndex].alignType = 1;
	    readList[validIndex].pos[1] = selfPos[0];
	    strcpy(readList[validIndex].seq[1],seq[0]);
	    rightAlignNum++;
	    validIndex++;
	 }
      }
      else if(strncmp(cigar[1],testCigar,index+1) == 0)
      {
	 initialReadStruct();
	 if(flag[1] == 153 || flag[1] == 147)
	 {
	    readList[validIndex].seq[1] = new char [readLength+1];
	    readList[validIndex].seq[0] = NULL;
	    readList[validIndex].alignType = 1;
	    readList[validIndex].pos[1] = selfPos[1];
	    strcpy(readList[validIndex].seq[1],seq[1]);
	    rightAlignNum++;
	    validIndex++;
	 }
	 else if(flag[1] == 137 || flag[1] == 163)
	 {
	    readList[validIndex].seq[0] = new char [readLength+1];
	    readList[validIndex].seq[1] = NULL;
	    readList[validIndex].alignType = 2;
	    readList[validIndex].pos[0] = selfPos[1];
	    strcpy(readList[validIndex].seq[0],seq[1]);
	    validIndex++;
	 }
      }
   }
}

void initialReadStruct()
{
   readList[validIndex].alignType = 0;
   readList[validIndex].pos[0] = 0;
   readList[validIndex].pos[1] = 0;
   readList[validIndex].mappingDist = 0;
   readList[validIndex].seq[0] = NULL;
   readList[validIndex].seq[1] = NULL;
}

void createSNPmatrix()
{
   int start = 0, blockCount = 0, startSNP[2], endSNP[2], snpNum, leftMost, rightMost, startSV;
   int i, j, k, l, index;
   char allele;
   bool left, right, sv_flag;

   snpMatrix = new int* [validIndex];
   for(i=rightAlignNum;i<validIndex;i++)
   {
      while(snpSite[start].getPos() < readList[i].pos[0] && start < snpCount) 
      {
	 if(snpSet[start] == -1) snpSet[start] = blockCount;
	 start++;
	 if(snpSet[start] == -1) blockCount++;
      }
      if(readList[i].alignType == 3)
      {
	 left = between(start,readList[i].pos[0],readList[i].pos[0]+readLength-1,&startSNP[0],&endSNP[0]);
	 right = between(start,readList[i].pos[1],readList[i].pos[1]+readLength-1,&startSNP[1],&endSNP[1]);
	 if(endSNP[1] <= endSNP[0]) right = false;
	 if(left && right)
	 {
	    if(startSNP[1] <= endSNP[0]) startSNP[1] = endSNP[0] + 1;
	    snpNum = (endSNP[0] - startSNP[0]) + (endSNP[1] - startSNP[1]) + 2;
	    snpMatrix[i] = new int [2*(snpNum+1)+1];
	    snpMatrix[i][0] = snpNum;
	    snpMatrix[i][1] = 0;
	    snpMatrix[i][2] = -1;
	    index = 3;
	    for(j=startSNP[0];j<=endSNP[1];j++) snpSet[j] = blockCount;
	    for(j=0;j<2;j++)
	    {
	       for(k=startSNP[j];k<=endSNP[j];k++)
	       {
		  if(snpSite[k].getAllele(0) == snpSite[k].getAllele(1)){
		     //cout << snpSite[k].getAllele(0) <<  snpSite[k].getAllele(1) << "\n" << endl ;
		     cout << snpSite[k].getPos()-readList[i].pos[j]+offset << endl;
		     allele = readList[i].seq[j][snpSite[k].getPos()-readList[i].pos[j]+offset];
		     if((snpSite[k].getPos()-readList[i].pos[j]+offset)>=120) break;
		     snpMatrix[i][index++] = allele;
		     snpMatrix[i][index++] = k;
		     cout << readList[i].pos[j] << snpSite[k].getAllele(0) << snpSite[k].getAllele(1) << allele << k << endl;
		     for(l=0;l<120;l++) {cout << readList[i].seq[j][l] ;}
		     cout << endl;
		  }else{
		     allele = readList[i].seq[j][snpSite[k].getPos()-readList[i].pos[j]];
		     snpMatrix[i][index++] = allele;
		     snpMatrix[i][index++] = k;
		  }
	       }
	    }
	 }
	 else if(left)
	 {
	    snpNum = (endSNP[0] - startSNP[0]) + 1;
	    snpMatrix[i] = new int [2*(snpNum+1)+1];
	    snpMatrix[i][0] = snpNum;
	    snpMatrix[i][1] = 0;
	    snpMatrix[i][2] = -1;
	    index = 3;
	    for(j=startSNP[0];j<=endSNP[0];j++) snpSet[j] = blockCount;
	    for(j=startSNP[0];j<=endSNP[0];j++)
	    {
	       allele = readList[i].seq[0][snpSite[j].getPos()-readList[i].pos[0]];
	       snpMatrix[i][index++] = allele;
	       snpMatrix[i][index++] = j;
	    }
	 }
	 else if(right)
	 {
	    snpNum = (endSNP[1] - startSNP[1]) + 1;
	    snpMatrix[i] = new int [2*(snpNum+1)+1];
	    snpMatrix[i][0] = snpNum;
	    snpMatrix[i][1] = 0;
	    snpMatrix[i][2] = -1;
	    index = 3;
	    for(j=startSNP[1];j<=endSNP[1];j++) snpSet[j] = blockCount;
	    for(j=startSNP[1];j<=endSNP[1];j++)
	    {
	       allele = readList[i].seq[1][snpSite[j].getPos()-readList[i].pos[1]];
	       snpMatrix[i][index++] = allele;
	       snpMatrix[i][index++] = j;
	    }
	 }
	 else
	 {
	    snpMatrix[i] = new int [1];
	    snpMatrix[i][0] = 0;
	 }
      }
      else if(readList[i].alignType == 2)
      {
	 left = between(start,readList[i].pos[0],readList[i].pos[0]+readLength-1,&startSNP[0],&endSNP[0]);
	 if(left)
	 {
	    sv_flag = false;
	    snpNum = (endSNP[0] - startSNP[0]) + 1;
	    snpMatrix[i] = new int [2*(snpNum+1)+1];
	    snpMatrix[i][0] = snpNum;
	    leftMost = readList[i].pos[0] + readLength;
	    rightMost = readList[i].pos[0] + insertSize + variance;
	    startSV = 0;
	    while(svSite[startSV].getPos(1) < leftMost && startSV < svCount) startSV++;
	    while(svSite[startSV].getPos(0) <= rightMost && startSV < svCount)
	    {
	       sv_flag = false;
	       if(svSite[startSV].overlap(leftMost,rightMost))
	       {
		  if(svSite[startSV].getMinSNP() == -1 && svSite[startSV].getMaxSNP() == -1)
		  {
		     svSite[startSV].setMinSNP(startSNP[0]);
		     svSite[startSV].setMaxSNP(endSNP[0]);
		  }
		  else if(startSNP[0] < svSite[startSV].getMinSNP()) svSite[startSV].setMinSNP(startSNP[0]);
		  else if(endSNP[0] > svSite[startSV].getMaxSNP()) svSite[startSV].setMaxSNP(endSNP[0]);
		  sv_flag = true;
		  break;
	       }
	       startSV++;
	    }
	    if(sv_flag && startSV < svCount)
	    {
	       snpMatrix[i][1] = svSite[startSV].getType();
	       snpMatrix[i][2] = startSV;
	    }	
	    else
	    {
	       snpMatrix[i][1] = 0;
	       snpMatrix[i][2] = -1;
	    }
	    index = 3;
	    for(j=startSNP[0];j<=endSNP[0];j++) snpSet[j] = blockCount;
	    for(j=startSNP[0];j<=endSNP[0];j++)
	    {
	       allele = readList[i].seq[0][snpSite[j].getPos()-readList[i].pos[0]];
	       snpMatrix[i][index++] = allele;
	       snpMatrix[i][index++] = j;
	    }
	 }
	 else
	 {
	    snpMatrix[i] = new int [1];
	    snpMatrix[i][0] = 0;
	 }
      }
      else
      {
	 snpMatrix[i] = new int [1];
	 snpMatrix[i][0] = 0;
      }
   }
   
   blockCount++;
   cout << rightAlignNum << "\t" << snpCount << "\t" << blockCount<<endl;
   while(start < snpCount)
   {
      if(snpSet[start] = -1)
      {
		 snpSet[start] = blockCount;
		 blockCount++;
      }
      start++;
   }
   start = 0;
   
   for(i=0;i<rightAlignNum;i++) 
   {
      while(snpSite[start].getPos() < readList[i].pos[1] && start < snpCount) start++;
      right = between(start,readList[i].pos[1],readList[i].pos[1]+readLength-1,&startSNP[1],&endSNP[1]);
      if(right)
      {
	 sv_flag = false;
	 snpNum = (endSNP[1] - startSNP[1]) + 1;
	 snpMatrix[i] = new int [2*(snpNum+1)+1];
	 snpMatrix[i][0] = snpNum;
	 leftMost = readList[i].pos[1] + readLength - insertSize - variance;
	 rightMost = readList[i].pos[1] - 1;
	 startSV = 0;
	 while(svSite[startSV].getPos(1) < leftMost && startSV < svCount) startSV++;
	 while(svSite[startSV].getPos(0) <= rightMost && startSV < svCount)
	 {
	    sv_flag = false;
	    if(svSite[startSV].overlap(leftMost,rightMost))
	    {
	       if(svSite[startSV].getMinSNP() == -1 && svSite[startSV].getMaxSNP() == -1)
	       {
		  svSite[startSV].setMinSNP(startSNP[1]);
		  svSite[startSV].setMaxSNP(endSNP[1]);
	       }
	       else if(startSNP[1] < svSite[startSV].getMinSNP()) svSite[startSV].setMinSNP(startSNP[1]);
	       else if(endSNP[1] > svSite[startSV].getMaxSNP()) svSite[startSV].setMaxSNP(endSNP[1]);
	       sv_flag = true;
	       break;
	    }
	    startSV++;
	 }
	 if(sv_flag && startSV < svCount)
	 {
	    snpMatrix[i][1] = svSite[startSV].getType();
	    snpMatrix[i][2] = startSV;
	 }
	 else
	 {
	    snpMatrix[i][1] = 0;
	    snpMatrix[i][2] = -1;
	 }
	 index = 3;
	 for(j=startSNP[1]+1;j<=endSNP[1];j++)
	 {
	    if(snpSet[j] != snpSet[startSNP[1]])
	    {
	       int tmp = snpSet[j];
	       for(k=j;snpSet[k]==tmp;k++) snpSet[k] = snpSet[startSNP[1]];
	    }
	 }
	 blockCount = snpSet[startSNP[1]];
	 for(j=startSNP[1];j<=endSNP[1];j++)
	 {
	    allele = readList[i].seq[1][snpSite[j].getPos()-readList[i].pos[1]];
	    snpMatrix[i][index++] = allele;
	    snpMatrix[i][index++] = j;
	 }
      }
      else
      {
	 snpMatrix[i] = new int [1];
	 snpMatrix[i][0] = 0;
      }
   }
}

bool between(int start, int leftMost, int rightMost, int *startSNP, int *endSNP)
{
   *startSNP = -1;
   *endSNP = -1;
   for(int i=start;i<snpCount;i++)
   {
      if(snpSite[i].getPos() >= leftMost && snpSite[i].getPos() <= rightMost)
      {
	 if(*startSNP == -1) *startSNP = i;
	 *endSNP = i;
      }
      if(snpSite[i].getPos() > rightMost) break;
   }
   if(*startSNP != -1) return true;
   else return false;
}

void printMatrix()
{
   int i, j, dice, cur;
   char *nucleotide = new char [4];
   FILE *fout;

   fout = fopen("tmpMatrix","w");
   
   for(i=0;i<validIndex;i++)
   {
      if(snpMatrix[i][0] == 0) continue;
      fprintf(fout,"%d %d\t",snpMatrix[i][1],snpMatrix[i][2]);

      for(j=0;j<snpMatrix[i][0];j++) 
      {
	 if(rand()%100 < errorRate)
	 {
	    if(snpMatrix[i][3+2*j] == 'A') cur = 0;
	    else if(snpMatrix[i][3+2*j] == 'C') cur = 1;
	    else if(snpMatrix[i][3+2*j] == 'G') cur = 2;
	    else if(snpMatrix[i][3+2*j] == 'T') cur = 3;
	    do
	    {
	       dice = rand()%4;
	    }
	    while(dice == cur);
	    if(dice == 0) snpMatrix[i][3+2*j] = 'A';
	    else if(dice == 1) snpMatrix[i][3+2*j] = 'C';
	    else if(dice == 2) snpMatrix[i][3+2*j] = 'G';
	    else if(dice == 3) snpMatrix[i][3+2*j] = 'T';
	 }
	 fprintf(fout,"%c %d ",snpMatrix[i][3+2*j],snpMatrix[i][4+2*j]);
      }
      fprintf(fout,"\n");
   }
   fclose(fout);
}

void printCluster()
{
   int i, j, curSet = snpSet[0], count = 0, min, max;
   FILE *fout;

   for(i=0;i<svCount;i++)
   {
      if(snpSet[svSite[i].getMinSNP()] != snpSet[svSite[i].getMaxSNP()])
      {
	 min = snpSet[svSite[i].getMinSNP()];
	 max = snpSet[svSite[i].getMaxSNP()];
	 for(j=svSite[i].getMinSNP()+1;j<snpCount;j++)
	 {
	    if(snpSet[j] > max) break;
	    if(snpSet[j] > min && snpSet[j] <= max) snpSet[j] = min;
	 }
      }
   }

   curSet = snpSet[0];
   fout = fopen("cluster.log","w");
   for(int i=0;i<snpCount;i++)
   {
      if(snpSet[i] != curSet)
      {
	 fprintf(fout,"\n%d ",i);
	 curSet = snpSet[i];
      }
      else fprintf(fout,"%d ",i);
   }
   fprintf(fout,"\n");
   fclose(fout);
}
