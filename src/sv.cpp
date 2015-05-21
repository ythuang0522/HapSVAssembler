#include <cstdlib>
#include <cstdio>
#include "SV.h"

SV::SV() 
{
   besideSnpList = new int [20];
   minSNP = -1;
   maxSNP = -1;
}

SV::~SV() 
{
   delete [] besideSnpList;
}

void SV::insertSV(int leftPos, int rightPos, int size, int type)
{
   position[0] = leftPos;
   position[1] = rightPos;
   svSize = size;
   svType = type;
}

int SV::getPos(int index)
{
   if(index != 0 && index != 1)
   {
      printf("ERROR PARAMETER! int SV::getPos(int index)\n");
      printf("Input index value: %d\n",index);
      printf("Structural variation exactly has 2 breakpoints\n\n");
      exit(0);
   }
   return position[index];
}

int SV::getSize()
{
   return svSize;
}

int SV::getType()
{
   return svType;
}

bool SV::overlap(int leftMost, int rightMost)
{
   if(rightMost < position[0]) return false;
   if(leftMost > position[1]) return false;
   return true;
}

void SV::printSVinfo()
{
   printf("%d\t%d\t[%d,%d]\n",getType(),svSize,position[0],position[1]);
}

void SV::insertBesideSNP(int snpIndex)
{
   besideSnpList[besideSnpNum] = snpIndex;
   besideSnpNum++;
}

int SV::getBesideSnpNum()
{
   return besideSnpNum;
}

int SV::getBesideSNP(int index)
{
   if(index < 0 || index >= besideSnpNum)
   {
      printf("ERROR PARAMETER! int SV::getBesideSNP(int index)\n");
      printf("Input index value: %d\n",index);
      printf("Only %d SNPs beside this SV\n\n",besideSnpNum);
      exit(0);
   }
   return besideSnpList[index];
}

void  SV::setMinSNP(int index)
{
   minSNP = index;
}

void SV::setMaxSNP(int index)
{
   maxSNP = index;
}

int SV::getMinSNP()
{
   return minSNP;
}

int SV::getMaxSNP()
{
   return maxSNP;
}
