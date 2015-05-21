#include <cstdlib>
#include <cstdio>
#include "SNP.h"

SNP::SNP() {}

SNP::~SNP() {}

void SNP::insertSNP(int pos, char allele1, char allele2)
{
   position = pos;
   allele[0] = allele1;
   allele[1] = allele2;
}

int SNP::getPos()
{
   return position;
}

char SNP::getAllele(int index)
{
   if(index != 0 && index != 1)
   {
      printf("ERROR PARAMRTER! char SNP::getAllele(int index)\n");
      printf("Input index value: %d\n",index);
      printf("Diploid genome has at most 2 alleles at a SNP site\n\n");
      exit(0);
   }
   return allele[index];
}

void SNP::printSNPinfo()
{
   printf("%d\t%c/%c\n",position,allele[0],allele[1]);
}
