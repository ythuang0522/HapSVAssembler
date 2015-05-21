#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <vector>
#include "BINDINGSITE.h"
#include "FRAGMENT.h"
#include "CLUSTER.h"
using namespace std;

class CHROMOSOME
{
   public:
      CHROMOSOME();
      ~CHROMOSOME();
      void insertInfo(FRAGMENT *srcFragListPtr, BINDINGSITE *srcBindingListPtr, int bind_num, int read_num, int* pheno, int* geno, int** srcSVList, int *sv_num);
      void initial();
      void geno2pheno();
      void pheno2geno();
      int maxElement(int *array);
      void evaluate();
      void constructHaplotype(int set, int index, int major);
      void setGene(int index, int allele);
      int getGene(int index);
      void flipGene(int index);
      int getChrSize();
      int getErrorCorrectionNum();
      char* getChr();
      char* getHaplotype(int index);
      void printErrorList();
      void updateSVHaplotype();
      int getErrorListSize();
      int getNormalListSize();
      int getErrorSite(int index);
      int getNormalSite(int index);
      bool checkChr();
      CHROMOSOME& operator=(const CHROMOSOME& obj);
   private:
      int *chr;
      int *chrGenotype;
      int *errorList;
      int *normalList;
      int *groupAllele[2];
      int *mapping2geno;
      int *mapping2pheno;
      char *chrString;
      char *haplotype[2];
      int **svList;
      BINDINGSITE *bind;
      FRAGMENT *fragListPtr;
      int *svNum;
      int bindNum;
      int snpNum;
      int readNum;
      int readNumAll;
      int errorReadNum;
      int normalReadNum;
      int errorCorrect;
};

#endif
