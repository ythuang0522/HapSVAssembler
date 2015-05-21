#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cstring>
using namespace std;

struct BP
{
   int pos;
   int side;
   bool used;
};

struct READ
{
   short int orientation;
   int leftPos;
   int rightPos;
   int mappingDist;
};

struct INDEL
{
   int type;
   int size;
   int avgFromMD;
   int leftBP;
   int rightBP;
   int count;
   bool bpRead;
};

struct INVERSION
{
   int size;
   int leftBP;
   int rightBP;
   int leftRegion[2];
   int rightRegion[2];
   int count;
   bool bpRead;
};

static int usage();
void readSamSummary();
void detectSV();
void detectInversion(int leftMost, int rightMost, bool side, int left, int right);
void detectIndel(int leftMost, int rightMost, int type, int svSize);
void readBP();
void preciseBP();
void outputSV();
int comparePos(const void *s1, const void *s2) {return ((*(struct READ*)s1).leftPos - (*(struct READ*)s2).leftPos);}
int compareIndelPos(const void *s1, const void *s2) {return ((*(struct INDEL*)s1).leftBP - (*(struct INDEL*)s2).leftBP);} 

int insertSize = 0;
int variance = 0;
int readLength = 0;
int row = 0;

int readNum = 0;
struct READ *readList;
int indelNum = 0;
struct INDEL *indelList;
int inversionNum = 0;
struct INVERSION *inversionList;
int bpNum = 0;
int totalBP = 0;
struct BP *bpList;
int curInversion = 0;
int curIndel = 0;
int normalOrientationNum = 0;
int svNum = 0;

void readSamSummary()
{
   ifstream fin;
   char buf[1024], cigar[2][256], *testCigar, digit[3];
   int flag[2], iSize[2], pos[2], mPos[2], tmpLength = readLength, i, index = 0, num;

   while(tmpLength/10 != 0)
   {
      digit[index] = tmpLength%10 + 48;
      tmpLength = tmpLength/10;
      index++;
   }
   digit[index] = tmpLength + 48;

   testCigar = new char [index+2];
   for(i=index;i>=0;i--) testCigar[index-i] = digit[i];
   testCigar[index+1] = 'M';
   testCigar[index+2] = '\0';

   fin.open("paired-align.dist");
   readList = new struct READ [row];

   //fix a bug of unpaired read in sam
   for(i=0;i<row;i+=2)
   {
      fin.getline(buf,1024);
      strcpy(cigar[0],strtok(buf,"\t"));
      flag[0]  = atoi(strtok(NULL,"\t"));
      iSize[0] = atoi(strtok(NULL,"\t"));
      pos[0]   = atoi(strtok(NULL,"\t")) - 1;
      mPos[0]  = atoi(strtok(NULL,"\n")) - 1;
      fin.getline(buf,1024);
      strcpy(cigar[1],strtok(buf,"\t"));
      if(strcmp(cigar[0],testCigar) != 0 || strcmp(cigar[1],testCigar) != 0) continue;
      flag[1]  = atoi(strtok(NULL,"\t"));
      iSize[1] = atoi(strtok(NULL,"\t"));
      pos[1]   = atoi(strtok(NULL,"\t")) - 1;
      mPos[1]  = atoi(strtok(NULL,"\n")) - 1;
      if((flag[0] == 99 && flag[1] == 147) || (flag[0] == 97 && flag[1] == 145)) 
      {
		 readList[readNum].orientation = 2; //(+,-)
		 normalOrientationNum++;
      }
      else if((flag[0] == 83 && flag[1] == 163) || (flag[0] == 81 && flag[1] == 161)) 
      {
		 readList[readNum].orientation = 2; //(-,+)
		 normalOrientationNum++;
      }
      else if(flag[0] == 65 && flag[1] == 129) readList[readNum].orientation = 3; //(+,+)
      else if(flag[0] == 113 && flag[1] == 177) readList[readNum].orientation = 0; //(-,-)
      else continue;
      if(iSize[0] > iSize[1])
      {
		 readList[readNum].leftPos = pos[0];
		 readList[readNum].rightPos = mPos[0];
		 readList[readNum].mappingDist = iSize[0];
      }
      else
      {
		 readList[readNum].leftPos = pos[1];
		 readList[readNum].rightPos = mPos[1];
		 readList[readNum].mappingDist = iSize[1];
      }
      readNum++;
   }
   fin.close();
   qsort(readList,readNum,sizeof(struct READ),comparePos);
}

void detectSV()
{
   int i, leftMost, rightMost, left, right, shift = readLength/10, type, size;
   bool side;
   inversionList = new struct INVERSION [readNum-normalOrientationNum];
   indelList = new struct INDEL [normalOrientationNum];

   for(i=0;i<readNum;i++)
   {
      if(readList[i].orientation == 3 || readList[i].orientation == 0) //discordant reads: match orientation factor
      {
	 if(readList[i].orientation == 3) //(+,+)
	 {
	    leftMost = readList[i].leftPos + readLength - shift;
	    rightMost = readList[i].leftPos + insertSize;
	    if(rightMost < readList[i].rightPos + readLength) rightMost = readList[i].rightPos + readLength;
	    side = false; //left breakpoint region
	 }
	 else //(-,-)
	 {
	    rightMost = readList[i].rightPos + shift - 1;
	    leftMost = readList[i].rightPos + readLength - insertSize;
	    if(leftMost > readList[i].leftPos) leftMost = readList[i].leftPos;
	    side = true; //right breakpoint region
	 }
	 left = readList[i].leftPos + readLength - shift;
	 right = readList[i].rightPos + shift - 1;
	 detectInversion(leftMost,rightMost,side,left,right);
      }
      else
      {
	 if(abs(readList[i].mappingDist - insertSize) >= 2*variance) //discordant reads: match distance factor
	 {
	    leftMost = readList[i].leftPos + readLength - shift;
	    rightMost = readList[i].rightPos + shift - 1;
	    (readList[i].mappingDist < insertSize)?(type = 1):(type = 0);
	    size = abs(readList[i].mappingDist - insertSize);
	    detectIndel(leftMost,rightMost,type,size);
	 }
      }
   }
}

void detectInversion(int leftMost, int rightMost, bool side, int left, int right)
{
   if(inversionNum == 0) //the first candidate discordant reads
   {
      inversionList[inversionNum].leftBP = leftMost;
      inversionList[inversionNum].rightBP = rightMost;
      inversionList[inversionNum].count = 1;
      inversionList[inversionNum].bpRead = false;
      inversionList[inversionNum].leftRegion[0] = 0;
      inversionList[inversionNum].leftRegion[1] = 0;
      inversionList[inversionNum].rightRegion[0] = 0;
      inversionList[inversionNum].rightRegion[1] = 0;
      if(!side) //left breakpoint region
      {
	 inversionList[inversionNum].leftRegion[0] = left;
	 inversionList[inversionNum].leftRegion[1] = right;
      }
      else //right breakpoint region
      {
	 inversionList[inversionNum].rightRegion[0] = left;
	 inversionList[inversionNum].rightRegion[1] = right;
      }
      inversionNum++;
   }
   else
   {
      while(leftMost > inversionList[curInversion].rightBP && curInversion < inversionNum) curInversion++;
      if(leftMost > inversionList[curInversion].rightBP || curInversion == inversionNum) //without overlapping
      {
	 inversionList[inversionNum].leftBP = leftMost;
	 inversionList[inversionNum].rightBP = rightMost;
	 inversionList[inversionNum].count = 1;
	 inversionList[inversionNum].bpRead = false;
	 inversionList[inversionNum].leftRegion[0] = 0;
	 inversionList[inversionNum].leftRegion[1] = 0;
	 inversionList[inversionNum].rightRegion[0] = 0;
	 inversionList[inversionNum].rightRegion[1] = 0;
	 if(!side) //left breakpoint region
	 {
	    inversionList[inversionNum].leftRegion[0] = left;
	    inversionList[inversionNum].leftRegion[1] = right;
	 }
	 else //right breakpoint region
	 {
	    inversionList[inversionNum].rightRegion[0] = left;
	    inversionList[inversionNum].rightRegion[1] = right;
	 }
	 inversionNum++;
      }
      else
      {
	 if(leftMost < inversionList[curInversion].leftBP) inversionList[curInversion].leftBP = leftMost;
	 if(rightMost > inversionList[curInversion].rightBP) inversionList[curInversion].rightBP = rightMost;
	 if(!side) //left breakpoint region
	 {
	    if(inversionList[curInversion].leftRegion[0] == 0 && inversionList[curInversion].leftRegion[1] == 0)
	    {
	       inversionList[curInversion].leftRegion[0] = left;
	       inversionList[curInversion].leftRegion[1] = left;
	    }
	    else
	    {
	       if(left > inversionList[curInversion].leftRegion[0]) inversionList[curInversion].leftRegion[0] = left;
	       if(right < inversionList[curInversion].leftRegion[1]) inversionList[curInversion].leftRegion[1] = right;
	    }
	 }
	 else
	 {
	    if(inversionList[curInversion].rightRegion[0] == 0 && inversionList[curInversion].rightRegion[1] == 0)
	    {
	       inversionList[curInversion].rightRegion[0] = left;
	       inversionList[curInversion].rightRegion[1] = right;
	    }
	    else
	    {
	       if(left > inversionList[curInversion].rightRegion[0]) inversionList[curInversion].rightRegion[0] = left;
	       if(right < inversionList[curInversion].rightRegion[1]) inversionList[curInversion].rightRegion[1] = right;
	    }
	 }
	 inversionList[curInversion].count++;
      }
   }
}

void detectIndel(int leftMost, int rightMost, int type, int svSize)
{
   int i;
   bool overlap = false;

   if(indelNum == 0)
   {
      indelList[indelNum].type = type;
      indelList[indelNum].leftBP = leftMost;
      indelList[indelNum].rightBP = rightMost;
      indelList[indelNum].count = 1;
      indelList[indelNum].avgFromMD = svSize;
      indelList[indelNum].bpRead = false;
      indelNum++;
   }
   else
   {
      for(i=0;i<indelNum;i++)
      {
	 if(leftMost > indelList[i].rightBP) continue;
	 if(rightMost < indelList[i].leftBP) break;
	 if(type == indelList[i].type /*&& (svSize>=(indelList[i].avgFromMD-variance) && svSize<=(indelList[i].avgFromMD+variance))*/)
	 {
	    if(leftMost > indelList[i].leftBP) indelList[i].leftBP = leftMost;
	    if(rightMost < indelList[i].rightBP) indelList[i].rightBP = rightMost;
	    indelList[i].count++;
	    indelList[i].avgFromMD = ((indelList[i].avgFromMD*(indelList[i].count-1))+svSize)/indelList[i].count;
	    overlap = true;
	    break;
	 }
      }
      if(!overlap)
      {
	 indelList[indelNum].type = type;
	 indelList[indelNum].leftBP = leftMost;
	 indelList[indelNum].rightBP = rightMost;
	 indelList[indelNum].count = 1;
	 indelList[indelNum].avgFromMD = svSize;
	 indelList[indelNum].bpRead = false;
	 indelNum++;
      }
      qsort(indelList,indelNum,sizeof(struct INDEL),compareIndelPos);
   }
}

void readBP()
{
   ifstream fin;
   char buf[1024], *tmp;

   fin.open("breakPoint.log");
   bpList = new struct BP [totalBP];

   while(fin.getline(buf,1024))
   {
      tmp = strtok(buf,"\t");
      bpList[bpNum].pos = atoi(tmp);
      tmp = strtok(NULL,"\n");
      bpList[bpNum].side = atoi(tmp);
      bpList[bpNum].used = false;
      if(bpNum > 0 && bpList[bpNum].pos >= bpList[bpNum-1].pos && bpList[bpNum].pos <= (bpList[bpNum-1].pos+variance))
      {
	 if(bpList[bpNum-1].side == 1 && bpList[bpNum].side == 2)
	 {
	    bpList[bpNum-1].side = 2;
	    bpList[bpNum].side = 1;
	 }
      }
      bpNum++;
   }
   fin.close();
}

void preciseBP()
{
   int i, bpIndex = 0, indelIndex = 0, curPos;

   for(i=0;i<inversionNum;i++)//dealing with inversion first
   {
      while(bpList[bpIndex].pos < inversionList[i].leftBP && bpIndex < bpNum) bpIndex++;
      while(bpList[bpIndex].pos >= inversionList[i].leftRegion[0] && bpList[bpIndex].pos <= inversionList[i].leftRegion[1])
      {
	 curPos = bpList[bpIndex].pos;
	 inversionList[i].leftRegion[0] = curPos;
	 inversionList[i].leftBP = curPos;
	 inversionList[i].bpRead = true;
	 bpList[bpIndex].used = true;
	 bpIndex++;
	 while(bpList[bpIndex].pos < (curPos+variance) && bpIndex < bpNum)
	 {
	    bpList[bpIndex].used = true;
	    bpIndex++;
	 }
	 if(bpIndex == bpNum) break;
      }
      while(bpList[bpIndex].pos >= inversionList[i].rightRegion[0] && bpList[bpIndex].pos <= inversionList[i].rightRegion[1])
      {
	 curPos = bpList[bpIndex].pos;
	 inversionList[i].rightRegion[1] = curPos;
	 inversionList[i].rightBP = curPos;
	 inversionList[i].bpRead = true;
	 bpList[bpIndex].used = true;
	 bpIndex++;
	 while(bpList[bpIndex].pos < (curPos+variance) && bpIndex < bpNum)
	 {
	    bpList[bpIndex].used = true;
	    bpIndex++;
	 }
	 if(bpIndex == bpNum) break;
      }
//      if(bpIndex == bpNum) break;
      inversionList[i].leftBP = inversionList[i].leftRegion[0];
      inversionList[i].rightBP = inversionList[i].rightRegion[1];
      inversionList[i].size = inversionList[i].rightBP - inversionList[i].leftBP + 1;
   }
   for(i=0;i<bpNum;i++)
   {
      if(bpList[i].used == false)
      {
	 while(bpList[i].pos > indelList[indelIndex].rightBP && indelIndex < indelNum) indelIndex++;
	 if(bpList[i].pos >= indelList[indelIndex].leftBP && bpList[i].pos <= indelList[indelIndex].rightBP)
	 {
	    if(bpList[i].side == 2) indelList[indelIndex].leftBP = bpList[i].pos; //leftBP
	    else indelList[indelIndex].rightBP = bpList[i].pos; //rightBP
	    indelList[indelIndex].size = indelList[indelIndex].rightBP - indelList[indelIndex].leftBP + 1;
	    indelList[indelIndex].bpRead = true;
	 }
      }
   }
}

int main(int argc, char *argv[])
{
   double startT, endT, overallTime = 0, sectionTime = 0;
   char c;

   while((c = getopt(argc,argv,"b:z:v:l:r:h")) >= 0)
   {
      switch(c)
      {
	 case 'b': totalBP = atoi(optarg); break;
	 case 'z': insertSize = atoi(optarg); break;
	 case 'v': variance = atoi(optarg); break;
	 case 'l': readLength = atoi(optarg); break;
	 case 'r': row = atoi(optarg); break;
	 case 'h': usage(); return 0;
	 default : usage(); return 1;
      }
   }
   if(optind != 11) return usage();

   startT = clock();
   readSamSummary();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[svDetector_readSamSummary] time elapses: %.2lf sec\n",sectionTime);
   cout << "[svDetector_readSamSummary] SAM has been loaded." << endl;

   startT = clock();
   detectSV();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[svDetector_detectSV] time elapses: %.2lf sec\n",sectionTime);
   cout << "[svDetector_detectSV] overall " << indelNum+inversionNum << " potiential SVs have been detected." << endl;

   startT = clock();
   readBP();
   preciseBP();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[svDetector_preciseBP] time elapses: %.2lf sec\n",sectionTime);
   cout << "[svDetector_preciseBP] precising the breakpoints of SVs." << endl;

   outputSV();

   printf( "[svDetector] overall time elapses: ");
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

void outputSV()
{
   FILE *fout;

   fout = fopen("Detected_SV_With_Count.log","w");

   for(int i=0;i<indelNum;i++)
   {
      if(indelList[i].count >= 5 || indelList[i].bpRead == true)
      {
	 if(indelList[i].type == 0) //deletion
	 {
	    indelList[i].size = indelList[i].rightBP - indelList[i].leftBP + 1;
	    if(indelList[i].size >= 5000) continue;
	    if(abs(indelList[i].avgFromMD - indelList[i].size) <= 2*variance && indelList[i].size >= 2*variance)
	    {
	       fprintf(fout,"Deletion\t");
	       fprintf(fout,"%d\t%d\t%d\tcount: %d\n",indelList[i].size,indelList[i].leftBP,indelList[i].rightBP,indelList[i].count);
	    }
	 }
	 else //insertion
	 {
	    indelList[i].size = indelList[i].avgFromMD;
	    if(indelList[i].size >= 5000) continue;
	    fprintf(fout,"Insertion\t");
	    fprintf(fout,"%d\t%d\t%d\tcount: %d\n",indelList[i].size,indelList[i].leftBP,indelList[i].rightBP,indelList[i].count);
	 }
	 svNum++;
      }
   }
   for(int i=0;i<inversionNum;i++)
   {
      if(inversionList[i].leftBP == 0 || inversionList[i].rightBP == 0) continue;
      if(inversionList[i].size < 0 || inversionList[i].size >= 5000) continue;
      if(inversionList[i].count >= 3)
      {
	 fprintf(fout,"Inversion\t%d\t%d\t%d\t",inversionList[i].size,inversionList[i].leftBP,inversionList[i].rightBP);
	 fprintf(fout,"count: %d\n",inversionList[i].count);
	 svNum++;
      }
   }
   fclose(fout);
}

static int usage()
{
   cout << "\nVersion 1.1: released on November 3rd, 2011" << endl;
   cout << "\nUsage: svDetector -z <int> -v <int> -l <int> -r <int> -b <int>" << endl;
   cout << "  -z <int> insert size of paired-end reads" << endl;
   cout << "  -v <int> standard deviation of insert size" << endl;
   cout << "  -l <int> number of known nucleotides on each ends" << endl;
   cout << "  -r <int> number of reads in SAM" << endl;
   cout << "  -b <int> number of detected breakpoints in breakpoint.log" << endl;
   cout << "  -h       help" << endl << endl;
   return 1;
}
