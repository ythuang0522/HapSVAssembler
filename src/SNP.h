#ifndef SNP_H
#define SNP_H

class SNP
{
   public:
      SNP();
      ~SNP();
      void insertSNP(int pos, char allele1, char allele2);
      int getPos();
      char getAllele(int index);
      void printSNPinfo();

   private:
      int position;
      char allele[2];
};

#endif
