#include "CLUSTER.h"
#include <iostream>
#include <cstdlib>
using namespace std;

CLUSTER::CLUSTER() 
{
   clusterSize = 0;
   fragNum = 0;
   fragStart = 0;
   pos = NULL;
}

CLUSTER::~CLUSTER() 
{
   delete [] pos;
}

void CLUSTER::setSize(int size)
{
   pos = new int [size];
   clusterSize = size;
}

void CLUSTER::setPos(int snpIndex, int snpPos)
{
   pos[snpIndex] = snpPos;
}

void CLUSTER::setFragNum(int num)
{
   fragNum = num;
}

void CLUSTER::setFragStart(int index)
{
   fragStart = index;
}

int CLUSTER::getPos(int index)
{
   if(index < 0 || index > clusterSize)
   {
      cout << "Index out of bound in CLUSTER::getPos(int index): index = " << index << endl;
      exit(0);
   }
   return pos[index];
}

int CLUSTER::getSize()
{
   return clusterSize;
}

int CLUSTER::getFragNum()
{
   return fragNum;
}

int CLUSTER::getFragStart()
{
   return fragStart;
}
