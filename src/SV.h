#ifndef SV_H
#define SV_H

class SV
{
   public:
      SV();
      ~SV();
      void insertSV(int leftPos, int rightPos, int size, int type);
      int getPos(int index);
      int getSize();
      int getType();
      bool overlap(int leftMost, int rightMost);
      void printSVinfo();
      void insertBesideSNP(int snpIndex);
      int getBesideSnpNum();
      int getBesideSNP(int index);
      int getMinSNP();
      void setMinSNP(int index);
      int getMaxSNP();
      void setMaxSNP(int index);

   private:
      int position[2];
      int svSize;
      int svType;
      int *besideSnpList;
      int besideSnpNum;
      int minSNP;
      int maxSNP;
};

#endif
