#include "BINDINGSITE.h"
#include <iostream>
#include <cstdlib>

BINDINGSITE::BINDINGSITE() 
{
   binding = NULL;
   value = 0;
   svNum = -1;
   bindingLength = 0;
}

BINDINGSITE::~BINDINGSITE() 
{
   delete [] binding;
}

void BINDINGSITE::allocate(int readNum)
{
   binding = new int [readNum];
}

void BINDINGSITE::insertSite(int index)
{
   binding[bindingLength] = index;
   bindingLength++;
}

int BINDINGSITE::getSize()
{
   return bindingLength;
}

int BINDINGSITE::getSite(int index)
{
   if(index < 0 || index >= bindingLength)
   {
      cout << "Index out of bound in int BINDINGSITE::getSite(int index): index = " << index << endl;
      exit(0);
   }
   return binding[index];
}

void BINDINGSITE::setValue(int val)
{
   value = val;
}

int BINDINGSITE::getValue()
{
   return value;
}

void BINDINGSITE::setSVNum(int num)
{
   svNum = num;
}

int BINDINGSITE::getSVNum()
{
   return svNum;
}
