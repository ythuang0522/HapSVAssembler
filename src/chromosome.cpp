#include "CHROMOSOME.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;

int maximum(int *src, char *pseudo, int index)
{
   int curMax = -1;
   int maxIndex = -1;

   for(int i=0;i<4;i++)
   {
      if(src[i] > curMax)
      {
	 curMax = src[i];
	 maxIndex = i;
      }
   }

   if(maxIndex == 0) pseudo[index] == 'A';
   else if(maxIndex == 1) pseudo[index] == 'C';
   else if(maxIndex == 2) pseudo[index] == 'G';
   else if(maxIndex == 3) pseudo[index] == 'T';
   else pseudo[index] == '-';
}

CHROMOSOME::CHROMOSOME()
{
   fragListPtr = NULL;
   svList = NULL;
   svNum = NULL;
   bind = NULL;
   chr = NULL;
   chrGenotype = NULL;
   haplotype[0] = NULL;
   haplotype[1] = NULL;
   errorList = NULL;
   normalList = NULL;
   chrString = NULL;
   bindNum = 0;
   snpNum = 0;
   readNum = 0;
   readNumAll = 0;
   normalReadNum = 0;
   errorReadNum = 0;
   errorCorrect = 0;
   groupAllele[0] = new int [4];
   groupAllele[1] = new int [4];
}

CHROMOSOME::~CHROMOSOME()
{
   delete [] chrGenotype;
   delete [] chr;
   delete [] chrString;
   delete [] errorList;
   delete [] normalList;
   delete [] haplotype[0];
   delete [] haplotype[1];
   delete [] groupAllele[0];
   delete [] groupAllele[1];
}

void CHROMOSOME::insertInfo(FRAGMENT *srcFragListPtr, BINDINGSITE *srcBindingListPtr, int bind_num, int read_num, int *pheno, int* geno, int** srcSVList, int *sv_num)
{
   fragListPtr = srcFragListPtr;
   svList = srcSVList;
   svNum = sv_num;
   bind = srcBindingListPtr;
   bindNum = bind_num;
   mapping2pheno = pheno;
   mapping2geno = geno;
   snpNum = fragListPtr[0].getSNPNum();
   readNumAll = read_num;
   chrGenotype = new int [readNumAll];
   chrString = new char [readNumAll+1];
   haplotype[0] = new char [snpNum+1];
   haplotype[1] = new char [snpNum+1];
   haplotype[0][snpNum] = '\0';
   haplotype[1][snpNum] = '\0';
   int bindRead = 0;
   for(int i=0;i<bindNum;i++) bindRead = bindRead + (bind[i].getSize() - 1);
   readNum = readNumAll - bindRead;
   chr = new int [readNum];
   for(int i=0;i<readNum;i++) chr[i] = -1;
   errorList = new int [readNum];
   normalList = new int [readNum];
}

void CHROMOSOME::initial()
{
   int i, j, k, set, readStart, snpIndex, index;
   int *noneDivide, noneDivideNum = 0;
   int set0, set1, none, tot;
   char allele;
   char *pseudoHaplotype[2];
   int **count;

   pseudoHaplotype[0] = new char [snpNum];
   pseudoHaplotype[1] = new char [snpNum];
   count = new int* [8];
   noneDivide = new int [readNumAll];
   for(i=0;i<8;i++) count[i] = new int [snpNum];

   for(i=0;i<8;i++) for(j=0;j<snpNum;j++) count[i][j] = 0;
   for(i=0;i<snpNum;i++) pseudoHaplotype[0][i] = '-';
   for(i=0;i<snpNum;i++) pseudoHaplotype[1][i] = '-';
   for(i=0;i<readNumAll;i++) chrGenotype[i] = -1;
   for(i=0;i<bindNum;i++)
   {
      set = rand()%2;
      bind[i].setValue(set);
      for(j=0;j<bind[i].getSize();j++)
      {
	 chrGenotype[bind[i].getSite(j)] = set;
	 for(k=0;k<snpNum;k++)
	 {
	    allele = fragListPtr[bind[i].getSite(j)].getAllele(k);
	    if(allele == '-') continue;
	    if(allele == 'A') count[2*set][k]++;
	    else if(allele == 'C') count[2*set+1][k]++;
	    else if(allele == 'G') count[2*set+2][k]++;
	    else if(allele == 'T') count[2*set+3][k]++;
	    if(pseudoHaplotype[set][k] == '-') pseudoHaplotype[set][k] = allele;
	    else
	    {
	       maximum(count[0],pseudoHaplotype[0],k);
	       maximum(count[4],pseudoHaplotype[1],k);
	    }
	 }
      }
   }
   if(readNumAll < 3) readStart = 0;
   else readStart = rand()%(readNumAll-2) + 1;
   for(i=0;i<snpNum;i++)
   {
      allele = fragListPtr[readStart].getAllele(i);
      if(allele != '-')
      {
	 snpIndex = i;
	 break;
      }
   }
   if(chrGenotype[readStart] == -1) 
   {
      set = rand()%2;
      chrGenotype[readStart] = set;
      for(i=snpIndex;i<snpNum;i++)
      {
	 allele = fragListPtr[readStart].getAllele(i);
	 if(allele != '-') 
	 {
	    if(allele == 'A') count[2*set][i]++;
	    else if(allele == 'C') count[2*set+1][i]++;
	    else if(allele == 'G') count[2*set+2][i]++;
	    else if(allele == 'T') count[2*set+3][i]++;
	    if(pseudoHaplotype[set][i] == '-') pseudoHaplotype[set][i] = allele;
	    else
	    {
	       maximum(count[0],pseudoHaplotype[0],i);
	       maximum(count[4],pseudoHaplotype[1],i);
	    }
	 }
      }
   }
   for(i=readStart+1;i<readNumAll;i++)
   {
      set0 = set1 = none = tot = 0;
      for(j=snpIndex;j<snpNum;j++)
      {
	 allele = fragListPtr[i].getAllele(j);
	 if(allele == '-') continue;
	 if(pseudoHaplotype[0][j] == '-' && pseudoHaplotype[1][j] == '-') none++;
	 else if(allele == pseudoHaplotype[0][j]) set0++;
	 else if(allele == pseudoHaplotype[1][j]) set1++;
	 else if(allele != pseudoHaplotype[0][j] && pseudoHaplotype[1][j] == '-') set1++;
	 else if(allele != pseudoHaplotype[1][j] && pseudoHaplotype[0][j] == '-') set0++;
	 else none++;
	 tot++;
      }
      if(tot == none)
      {
	 noneDivide[noneDivideNum] = i;
	 noneDivideNum++;
	 continue;
      }
      if(set0 > set1)
      {
	 for(j=snpIndex;j<snpNum;j++)
	 {
	    allele = fragListPtr[i].getAllele(j);
	    if(allele != '-')
	    {
	       if(allele == 'A') count[0][j]++;
	       else if(allele == 'C') count[1][j]++;
	       else if(allele == 'G') count[2][j]++;
	       else if(allele == 'T') count[3][j]++;
	       if(pseudoHaplotype[0][j] == '-') pseudoHaplotype[0][j] = allele;
	       else
	       {
		  maximum(count[0],pseudoHaplotype[0],j);
		  maximum(count[4],pseudoHaplotype[1],j);
	       }
	    }
	 }
	 chrGenotype[i] = 0;
      }
      else if(set1 > set0)
      {
	 for(j=snpIndex;j<snpNum;j++)
	 {
	    allele = fragListPtr[i].getAllele(j);
	    if(allele != '-')
	    {
	       if(allele == 'A') count[4][j]++;
	       else if(allele == 'C') count[5][j]++;
	       else if(allele == 'G') count[6][j]++;
	       else if(allele == 'T') count[7][j]++;
	       if(pseudoHaplotype[1][j] == '-') pseudoHaplotype[1][j] = allele;
	       else
	       {
		  maximum(count[0],pseudoHaplotype[0],j);
		  maximum(count[4],pseudoHaplotype[1],j);
	       }
	    }
	 }
	 chrGenotype[i] = 1;
      }
      else
      {
	 noneDivide[noneDivideNum] = i;
	 noneDivideNum++;
      }
   }
   for(i=readStart-1;i>=0;i--)
   {
      set0 = set1 = none = tot = 0;
      for(j=0;j<snpNum;j++)
      {
	 allele = fragListPtr[i].getAllele(j);
	 if(allele == '-') continue;
	 if(pseudoHaplotype[0][j] == '-' && pseudoHaplotype[1][j] == '-') none++;
	 else if(allele == pseudoHaplotype[0][j]) set0++;
	 else if(allele == pseudoHaplotype[1][j]) set1++;
	 else if(allele != pseudoHaplotype[0][j] && pseudoHaplotype[1][j] == '-') set1++;
	 else if(allele != pseudoHaplotype[1][j] && pseudoHaplotype[0][j] == '-') set0++;
	 else none++;
	 tot++;
      }
      if(tot == none)
      {
	 noneDivide[noneDivideNum] = i;
	 noneDivideNum++;
	 continue;
      }
      if(set0 > set1)
      {
	 for(j=0;j<snpNum;j++)
	 {
	    allele = fragListPtr[i].getAllele(j);
	    if(allele != '-')
	    {
	       if(allele == 'A') count[0][j]++;
	       else if(allele == 'C') count[1][j]++;
	       else if(allele == 'G') count[2][j]++;
	       else if(allele == 'T') count[3][j]++;
	       if(pseudoHaplotype[0][j] == '-') pseudoHaplotype[0][j] = allele;
	       else
	       {
		  maximum(count[0],pseudoHaplotype[0],j);
		  maximum(count[4],pseudoHaplotype[1],j);
	       }
	    }
	 }
	 chrGenotype[i] = 0;
      }
      else if(set1 > set0)
      {
	 for(j=0;j<snpNum;j++)
	 {
	    allele = fragListPtr[i].getAllele(j);
	    if(allele != '-')
	    {
	       if(allele == 'A') count[4][j]++;
	       else if(allele == 'C') count[5][j]++;
	       else if(allele == 'G') count[6][j]++;
	       else if(allele == 'T') count[7][j]++;
	       if(pseudoHaplotype[1][j] == '-') pseudoHaplotype[1][j] = allele;
	       else
	       {
		  maximum(count[0],pseudoHaplotype[0],j);
		  maximum(count[4],pseudoHaplotype[1],j);
	       }
	    }
	 }
	 chrGenotype[i] = 1;
      }
      else
      {
	 noneDivide[noneDivideNum] = i;
	 noneDivideNum++;
      }
   }
   for(i=0;i<noneDivideNum;i++)
   {
      index = noneDivide[i];
      set0 = set1 = none = tot = 0;
      for(j=0;j<snpNum;j++)
      {
	 allele = fragListPtr[index].getAllele(j);
	 if(allele == '-') continue;
	 if(pseudoHaplotype[0][j] == '-' && pseudoHaplotype[1][j] == '-') none++;
	 else if(allele == pseudoHaplotype[0][j]) set0++;
	 else if(allele == pseudoHaplotype[1][j]) set1++;
	 else if(allele != pseudoHaplotype[0][j] && pseudoHaplotype[1][j] == '-') set1++;
	 else if(allele != pseudoHaplotype[1][j] && pseudoHaplotype[0][j] == '-') set0++;
	 else none++;
	 tot++;
      }
      if(set0 > set1)
      {
	 for(j=0;j<snpNum;j++)
	 {
	    allele = fragListPtr[index].getAllele(j);
	    if(allele != '-')
	    {
	       if(allele == 'A') count[0][j]++;
	       else if(allele == 'C') count[1][j]++;
	       else if(allele == 'G') count[2][j]++;
	       else if(allele == 'T') count[3][j]++;
	       if(pseudoHaplotype[0][j] == '-') pseudoHaplotype[0][j] = allele;
	       else
	       {
		  maximum(count[0],pseudoHaplotype[0],j);
		  maximum(count[4],pseudoHaplotype[1],j);
	       }
	    }
	 }
	 chrGenotype[index] = 0;
      }
      else if(set1 > set0)
      {
	 for(j=0;j<snpNum;j++)
	 {
	    allele = fragListPtr[index].getAllele(j);
	    if(allele != '-')
	    {
	       if(allele == 'A') count[4][j]++;
	       else if(allele == 'C') count[5][j]++;
	       else if(allele == 'G') count[6][j]++;
	       else if(allele == 'T') count[7][j]++;
	       if(pseudoHaplotype[1][j] == '-') pseudoHaplotype[1][j] = allele;
	       else
	       {
		  maximum(count[0],pseudoHaplotype[0],j);
		  maximum(count[4],pseudoHaplotype[1],j);
	       }
	    }
	 }
	 chrGenotype[index] = 1;
      }
      else
      {
	 set = rand()%2;
	 chrGenotype[index] = set;
	 for(j=0;j<snpNum;j++)
	 {
	    allele = fragListPtr[index].getAllele(j);
	    if(pseudoHaplotype[chrGenotype[index]][j] == '-' && allele != '-') pseudoHaplotype[chrGenotype[index]][j] = allele;
	 }

      }
   }

   geno2pheno();

   delete [] pseudoHaplotype[0];
   delete [] pseudoHaplotype[1];
   for(i=0;i<4;i++) delete [] count[i];
   delete [] noneDivide;
}

void CHROMOSOME::geno2pheno()
{
   int index = bindNum;

   for(int i=0;i<bindNum;i++)
   {
      chr[i] = bind[i].getValue();
      for(int j=0;j<bind[i].getSize();j++) chrGenotype[bind[i].getSite(j)] = -1;
   }
   for(int i=0;i<readNumAll;i++)
   {
      if(chrGenotype[i] != -1)
      {
	 chr[index] = chrGenotype[i];
	 index++;
      }
   }
}

void CHROMOSOME::pheno2geno()
{
   int index = 0, set;

   for(int i=0;i<readNumAll;i++) chrGenotype[i] = -1;
   for(int i=0;i<bindNum;i++)
   {
      set = bind[i].getValue();
      for(int j=0;j<bind[i].getSize();j++) chrGenotype[bind[i].getSite(j)] = set;
   }
   for(int i=bindNum;i<readNum;i++)
   {
      while(chrGenotype[index] != -1) index++;
      chrGenotype[index] = chr[i];
      index++;
   }
}

int CHROMOSOME::maxElement(int *array)
{
   int max = 0, index = 0;

   for(int i=0;i<4;i++)
   {
      if(array[i] > max)
      {
	 max = array[i];
	 index = i;
      }
   }
   return index;
}

void CHROMOSOME::evaluate()
{
   pheno2geno();
   bool *error = new bool [readNum];
   for(int i=0;i<readNum;i++) error[i] = false;

   int i, j;
   char allele;

   for(i=0;i<snpNum;i++)
   {
      for(j=0;j<4;j++)
      {
	 groupAllele[0][j] = 0;
	 groupAllele[1][j] = 0;
      }
      for(j=0;j<readNumAll;j++)
      {
	 allele = fragListPtr[j].getAllele(i);
	 if(allele == '-') continue;
	 if(allele == 'A') groupAllele[chrGenotype[j]][0]++;
	 else if(allele == 'C') groupAllele[chrGenotype[j]][1]++;
	 else if(allele == 'G') groupAllele[chrGenotype[j]][2]++;
	 else if(allele == 'T') groupAllele[chrGenotype[j]][3]++;
      }
      int major[2];
      major[0] = maxElement(groupAllele[0]);
      major[1] = maxElement(groupAllele[1]);
      if(groupAllele[0][major[0]] == 0 && groupAllele[1][major[1]] == 0)
      {
	 constructHaplotype(0,i,5);
	 constructHaplotype(1,i,5);
      }
      else
      {
	 if(groupAllele[0][major[0]] == 0) major[0] = major[1]; 
	 if(groupAllele[1][major[1]] == 0) major[1] = major[0];
	 if(major[0] == major[1])
	 {
	    if(groupAllele[0][major[0]] >= groupAllele[1][major[1]])
	    {
	       groupAllele[1][major[1]] = -1;
	       major[1] = maxElement(groupAllele[1]);
	       if(groupAllele[1][major[1]] == 0) major[1] = major[0];
	    }
	    else if(groupAllele[0][major[0]] < groupAllele[1][major[1]])
	    {
	       groupAllele[0][major[0]] = -1;
	       major[0] = maxElement(groupAllele[0]);
	       if(groupAllele[0][major[0]] == 0) major[0] = major[1];
	    }
	 }
	 constructHaplotype(0,i,major[0]);
	 constructHaplotype(1,i,major[1]);
      }
   }
   errorCorrect = 0;
   errorReadNum = 0;
   normalReadNum = 0;
   for(i=0;i<readNumAll;i++)
   {
      bool errorFlag = false;
      for(j=0;j<snpNum;j++)
      {
	 allele = fragListPtr[i].getAllele(j);
	 if(allele == '-') continue;
	 if(haplotype[chrGenotype[i]][j] != allele)
	 {
	    errorCorrect++;
	    errorFlag = true;
	 }
      }
      if(errorFlag == true) error[mapping2pheno[i]] = true;
   }
   for(i=0;i<readNum;i++)
   {
      if(error[i] == true) 
      {
	 errorList[errorReadNum] = i;
	 errorReadNum++;
      }
      else
      {
	 normalList[normalReadNum] = i;
	 normalReadNum++;
      }
   }
   delete [] error;
}

void CHROMOSOME::constructHaplotype(int set, int index, int major)
{
   if(major == 0) haplotype[set][index] = 'A';
   else if(major == 1) haplotype[set][index] = 'C';
   else if(major == 2) haplotype[set][index] = 'G';
   else if(major == 3) haplotype[set][index] = 'T';
   else haplotype[set][index] = '-';
}

void CHROMOSOME::updateSVHaplotype()
{
   bool flag = false;

   if(bindNum != 0)
   {
      for(int i=0;i<bindNum;i++)
      {
//	 flag = false;

//	 for(int j=0;j<*svNum;j++)
//	 {
//	    if(svList[j][0] == bind[i].getSVNum()) 
//	    {
//	       svList[j][1] = bind[i].getValue();
//	       flag = true;
//	       break;
//	    }
//	 }
//	 if(flag == false)
//	 {	 
	    svList[*svNum][0] = bind[i].getSVNum();
	    svList[*svNum][1] = chr[i];
	    (*svNum)++;
//	 }
      }
   }
}

void CHROMOSOME::setGene(int index, int allele)
{
   if(allele != 0 && allele != 1)
   {
      cout << "ERROR! void CHROMOSOME::setGene(int index, int allele) in chromosome.cpp\nallele = " << allele << endl;
      cout << "allele domain = [0,1]" << endl;
      exit(0);
   }
   if(index < 0 || index > readNum)
   {
      cout << "ERROR! void CHROMOSOME::setGene(int index, int allele) in chromosome.cpp\nindex = " << index << endl;
      cout << "index domain = [0," << readNum << "]" << endl;
      exit(0);
   }
   chr[index] = allele;
}

int CHROMOSOME::getGene(int index)
{
   if(index < 0 || index > readNum)
   {
      cout << "ERROR! void CHROMOSOME::getGene(int index) in chromosome.cpp\nindex = " << index << endl;
      cout << "index domain = [0," << readNum << "]" << endl;
      exit(0);
   }
   return chr[index];
}

int CHROMOSOME::getChrSize()
{
   return readNum;
}

int CHROMOSOME::getErrorCorrectionNum()
{
   return errorCorrect;
}

char* CHROMOSOME::getChr()
{
   for(int i=0;i<readNum;i++) chrString[i] = chr[i] + 48;
   chrString[readNum] = '\0';
   return chrString;
}

char* CHROMOSOME::getHaplotype(int index)
{
   if(index != 0 && index != 1)
   {
      cout << "ERROR! void CHROMOSOME::getHaplotype(int index) in chromosome.cpp\nindex = " << index << endl;
      cout << "index domain = 0 or 1" << endl;
      exit(0);
   }
   return haplotype[index];
}

void CHROMOSOME::flipGene(int index)
{
   if(index < 0 || index > readNum)
   {
      cout << "ERROR! void CHROMOSOME::flipGene(int index) in chromosome.cpp\nindex = " << index << endl;
      cout << "index domain = [0," << readNum << "]" << endl;
      exit(0);
   }
   chr[index] = (chr[index]+1)%2;
}

void CHROMOSOME::printErrorList()
{
   if(errorReadNum == 0) return;
   for(int i=0;i<errorReadNum;i++) printf("%4d ",errorList[i]);
   printf("\n");
}

int CHROMOSOME::getErrorListSize()
{
   return errorReadNum;
}

int CHROMOSOME::getNormalListSize()
{
   return normalReadNum;
}

int CHROMOSOME::getErrorSite(int index)
{
   if(index < 0 || index >= errorReadNum)
   {
      cout << "ERROR! void CHROMOSOME::getErrorSite(int index)(int index) in chromosome.cpp\nindex = " << index << endl;
      cout << "index domain = [0," << errorReadNum << "]" << endl;
      exit(0);
   }
   return errorList[index];
}

int CHROMOSOME::getNormalSite(int index)
{
   if(index < 0 || index >= normalReadNum)
   {
      cout << "ERROR! void CHROMOSOME::getNormalSite(int index)(int index) in chromosome.cpp\nindex = " << index << endl;
      cout << "index domain = [0," << normalReadNum << "]" << endl;
      exit(0);
   }
   return normalList[index];
}

CHROMOSOME& CHROMOSOME::operator=(const CHROMOSOME& obj)
{
   if(this != &obj)
   {
      bind = obj.bind;
      mapping2pheno = obj.mapping2pheno;
      mapping2geno = obj.mapping2geno;
      fragListPtr = obj.fragListPtr;
      svList = obj.svList;
      svNum = obj.svNum;
      bindNum = obj.bindNum;
      snpNum = obj.snpNum;
      readNum = obj.readNum;
      readNumAll = obj.readNumAll;
      normalReadNum = obj.normalReadNum;
      errorReadNum = obj.errorReadNum;
      errorCorrect = obj.errorCorrect;

      delete [] chr;
      delete [] chrGenotype;
      delete [] errorList;
      delete [] normalList;
      delete [] chrString;
      delete [] groupAllele[0];
      delete [] groupAllele[1];
      delete [] haplotype[0];
      delete [] haplotype[1];

      chr = new int [readNum];
      chrGenotype = new int [readNumAll];
      errorList = new int [readNum];
      normalList = new int [readNum];
      chrString = new char [readNumAll+1];
      groupAllele[0] = new int [4];
      groupAllele[1] = new int [4];
      haplotype[0] = new char [snpNum+1];
      haplotype[1] = new char [snpNum+1];

      memcpy(chr,obj.chr,sizeof(int)*readNum);
      memcpy(errorList,obj.errorList,sizeof(int)*readNum);
      memcpy(normalList,obj.normalList,sizeof(int)*readNum);
      memcpy(haplotype[0],obj.haplotype[0],sizeof(char)*(snpNum+1));
      memcpy(haplotype[1],obj.haplotype[1],sizeof(char)*(snpNum+1));
   }
   return *this;
}

bool CHROMOSOME::checkChr()
{
   for(int i=0;i<readNum;i++)
   {
      if(chr[i] != 0 && chr[i] != 1) return false;
   }
   return true;
}
