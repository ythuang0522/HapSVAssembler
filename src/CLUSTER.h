#ifndef CLUSTER_H
#define CLUSTER_H

using namespace std;

class CLUSTER
{
   public:
      CLUSTER();
      ~CLUSTER();
      void setSize(int size);
      void setPos(int snpIndex, int snpPos);
      void setFragNum(int num);
      void setFragStart(int index);
      int getPos(int index);
      int getSize();
      int getFragNum();
      int getFragStart();
   private:
      int *pos;
      int clusterSize;
      int fragNum;
      int fragStart;
};

#endif
