#ifndef BINDINGSITE_H
#define BINDINGSITE_H

using namespace std;

class BINDINGSITE
{
   public:
      BINDINGSITE();
      ~BINDINGSITE();
      void allocate(int readNum);
      void insertSite(int index);
      int getSize();
      int getSite(int index);
      void setValue(int val);
      int getValue();
      void setSVNum(int num);
      int getSVNum();

   private:
      int value;
      int *binding;
      int svNum;
      int bindingLength;
};

#endif
