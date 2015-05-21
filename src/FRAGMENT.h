#ifndef FRAGMENT_H
#define FRAGMENT_H

#include <string>
using namespace std;

class FRAGMENT
{
   public:
      FRAGMENT();
      ~FRAGMENT();
      void insertContent(char *srcRead, int type, int num, int len);
      char getAllele(int index);
      char* getRead();
      int getSVType();
      int getSVNum();
      int getSNPNum();

   private:
      char *read;
      int svType;
      int svNum;
      int snpNum;
};

#endif
