#include "FRAGMENT.h"
#include <cstdlib>
#include <cstring>

FRAGMENT::FRAGMENT() 
{
   read = NULL;
   svType = -1;
   svNum = -1;
   snpNum = 0;
}

FRAGMENT::~FRAGMENT() 
{
   delete [] read;
}

void FRAGMENT::insertContent(char *srcRead, int type, int num, int len)
{
   snpNum = len;
   read = new char [snpNum+1];
   strncpy(read,srcRead,snpNum);
   read[snpNum] = '\0';
   svType = type;
   svNum = num;
}

char FRAGMENT::getAllele(int index)
{
   return read[index];
}

char* FRAGMENT::getRead()
{
   return read;
}

int FRAGMENT::getSVType()
{
   return svType;
}

int FRAGMENT::getSVNum()
{
   return svNum;
}

int FRAGMENT::getSNPNum()
{
   return snpNum;
}
