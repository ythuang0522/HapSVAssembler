#include <unistd.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <pthread.h>

#include "CLUSTER.h"
#include "FRAGMENT.h"
#include "CHROMOSOME.h"
#include "BINDINGSITE.h"
using namespace std;

CLUSTER *cluster;
FRAGMENT *fragment;
BINDINGSITE *bind;
CHROMOSOME *population;
CHROMOSOME *children;
pthread_t *thrd;
int **svList;
int svCount = 0;
int svNum = 0;
int *mapping2geno;
int *mapping2pheno;
int popSize = 10;
int generation = 5;
int mutateRate = 80;
int readCount = 0;
int snpCount = 0;
int clusterCount = 0;
int clusterIndex = 0;
int bindingCount = 0;
char *haplotype[2];

void swapChromosome(CHROMOSOME &variable1, CHROMOSOME &variable2)
{
   CHROMOSOME temp;
   temp = variable1;
   variable1 = variable2;
   variable2 = temp;
}
static int usage();
void readCluster();
void readFragment();
void checkBinding(int clusterNum);
void initialMapping(int clusterNum);
void pipeline();
void parentSelection(int *parent);
void crossOver(int *parent, int index);
void mutation(CHROMOSOME* pop);
void sortPop();
void sortChild();
void survival();
int compareErrorNum(const void *s1, const void *s2);
void* parallel(void* data);
void* parallel_init(void* data);
void printHaplotype();

int main(int argc, char *argv[])
{
   double startT, endT, overallTime = 0, sectionTime = 0;
   time_t start_tm, finish_tm; 
   char c;

   while((c = getopt(argc,argv,"s:r:c:V:p:g:m:h")) >= 0)
   {
      switch(c)
      {
	 case 's': snpCount = atoi(optarg); break;
	 case 'r': readCount = atoi(optarg); break;
	 case 'c': clusterCount = atoi(optarg); break;
	 case 'V': svCount = atoi(optarg); break;
	 case 'p': popSize = atoi(optarg); break;
	 case 'g': generation = atoi(optarg); break;
	 case 'm': mutateRate = atoi(optarg); break;
	 case 'h': usage(); return 0;
	 default : cout << "Unrecognized option: -" << c << endl; return 1;
      }
   }
   if(optind < 9) return usage();
   if(mutateRate < 0 || mutateRate > 100)
   {
      cout << "Inappropriate parameter -m " << mutateRate << endl;
      return usage();
   }

   srand(time((unsigned)0));
   printf( "[hapAssembler_arguments] %s -s %d -r %d -c %d ",argv[0],snpCount,readCount,clusterCount);
   printf( "-V %d -p %d -g %d -m %d\n",svCount,popSize,generation,mutateRate);

   startT = clock();
   cluster = new CLUSTER[clusterCount];
   readCluster();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[hapAssembler_readCluster] time elapses: %.2lf sec\n",sectionTime);
   cout << "[hapAssembler_readCluster] " << clusterCount << " clusters have been loaded." << endl;

   startT = clock();
   fragment = new FRAGMENT[readCount];
   readFragment();
   endT = clock();
   sectionTime = (endT-startT)/CLOCKS_PER_SEC;
   overallTime = overallTime + sectionTime;
   printf( "[hapAssembler_readFragment] time elapses: %.2lf sec\n",sectionTime);
   cout << "[hapAssembler_readFragment] " << readCount << " SNP fragments have been loaded." << endl;

   time(&start_tm); 
   pipeline();
   time(&finish_tm); 
   sectionTime = difftime(finish_tm,start_tm); 
   overallTime = overallTime + sectionTime;
   printf( "[hapAssembler_pipeline] time elapses: %.2lf sec\n",sectionTime);
   cout << "[hapAssembler_pipeline] GA-pipeline has been processed." << endl;

   printf( "[hapAssembler] overall time elapses: ");
   if(overallTime/3600 > 1)
   {
      printf( "%.0lf Hours ",overallTime/3600);
      overallTime = overallTime - floor(overallTime/3600)*3600;
   }
   if(overallTime/60 > 1)
   {
      printf( "%.0lf Minutes ",overallTime/60);
      overallTime = overallTime - floor(overallTime/60)*60;
   }
   printf( "%.2lf Seconds\n",overallTime);

   delete [] cluster;
   delete [] fragment;

   return 0;
}

static int usage()
{
   cout << "\nVersion 1.1: released on November 3rd, 2011" << endl;
   cout << "\nUsage: hapAssembler -s <int> -r <int> -c <int> -V <int> [option]" << endl;
   cout << "  -s <int> number of SNPs" << endl;
   cout << "  -r <int> number of rows in SNP matrix" << endl;
   cout << "  -V <int> number of structural variations" << endl;
   cout << "  -c <int> number of SV clusters" << endl << endl;
   cout << "  [option]" << endl;
   cout << "  -p <int> population size [Default 10]" << endl;
   cout << "  -g <int> generation for GA-pipeline [Default 5]" << endl;
   cout << "  -m <int> mutation rate(m/100) [Default 80]" << endl;
   cout << "  -h       help" << endl << endl;
   return 1;
}

void readCluster()
{
   ifstream fin;
   char buf[1000000];
   char *tmp;
   int *tmpList = new int [snpCount];
   int count; //number of SNP sites per cluster

   fin.open("cluster.log");
   for(int i=0;i<clusterCount;i++)
   {
      count = 0;
      fin.getline(buf,1000000);
      tmp = strtok(buf," \n");
      tmpList[count] = atoi(tmp);
      count++;
      while((tmp = strtok(NULL," \n")) != NULL)
      {
	 tmpList[count] = atoi(tmp);
	 count++;
      }
      cluster[i].setSize(count); //dynamically resize the pos list in class CLUSTER
      for(int j=0;j<count;j++) cluster[i].setPos(j,tmpList[j]);
   }
   fin.close();
   delete [] tmpList;
}

void readFragment()
{
   ifstream fin;
   char buf[5000], seq[5000];
   char *tmp, allele;
   int type, num, firstPos, len, i, j;
   int boundary[2], count = 0;
   bool flag = false;

   fin.open("snp_matrix.log");
   len = cluster[clusterIndex].getSize();
   boundary[0] = cluster[clusterIndex].getPos(0);
   boundary[1] = cluster[clusterIndex].getPos(len-1);
   cluster[clusterIndex].setFragStart(0);
   for(i=0;i<readCount;i++)
   {
      for(j=0;j<5000;j++) seq[j] = '-';
      fin.getline(buf,5000);
      tmp = strtok(buf," "); //SV type
      type = atoi(tmp);
      tmp = strtok(NULL," \t"); //SV number
      num = atoi(tmp);
      tmp = strtok(NULL," "); //first allele
      allele = *tmp;
      tmp = strtok(NULL," \n"); //the position of first allele
      firstPos = atoi(tmp);
      
      if(firstPos >= boundary[0] && firstPos <= boundary[1]) count++;
      else
      {
	 if(count != 0) cluster[clusterIndex].setFragNum(count);
	 else
	 {
	    cluster[clusterIndex].setFragNum(0);
	    cluster[clusterIndex].setFragStart(-1);
	 }
	 clusterIndex++;
	 cluster[clusterIndex].setFragStart(i);
	 len = cluster[clusterIndex].getSize();
	 boundary[0] = cluster[clusterIndex].getPos(0);
	 boundary[1] = cluster[clusterIndex].getPos(len-1);
	 if(firstPos >= boundary[0] && firstPos <= boundary[1]) count = 1;
	 else count = 0;
      }
      seq[firstPos-boundary[0]] = allele;
      while((tmp = strtok(NULL," \n")) != NULL)
      {
	 allele = *tmp;
	 tmp = strtok(NULL," \n");
	 seq[atoi(tmp)-boundary[0]] = allele;
      }
      fragment[i].insertContent(seq,type,num,len);
   }
   if(count != 0) cluster[clusterIndex].setFragNum(count);
   else
   {
      cluster[clusterIndex].setFragNum(0);
      cluster[clusterIndex].setFragStart(-1);
   }
   clusterIndex++;
   while(clusterIndex < clusterCount)
   {
      cluster[clusterIndex].setFragNum(0);
      cluster[clusterIndex].setFragStart(-1);
      clusterIndex++;
   }
   fin.close();
}

void checkBinding(int clusterNum)
{
   int start = cluster[clusterNum].getFragStart();
   int num = cluster[clusterNum].getFragNum();
   int end = start + num;
   bool flag = false;

   bind = NULL;
   bindingCount = 0;
   for(int i=start;i<end;i++)
   {
      if(fragment[i].getSVNum() != -1)
      {
	 if(bind == NULL)
	 {
	    bind = new BINDINGSITE [1000];
	    bind[0].allocate(num);
	    bind[0].setSVNum(fragment[i].getSVNum());
	    bind[0].insertSite(i-start);
	    bindingCount++;
	 }
	 else
	 {
	    flag = false;
	    for(int j=0;j<bindingCount;j++)
	    {
	       if(fragment[i].getSVNum() == bind[j].getSVNum())
	       {
		  bind[j].insertSite(i-start);
		  flag = true;
		  break;
	       }
	    }
	    if(flag == false)
	    {
	       bind[bindingCount].allocate(num);
	       bind[bindingCount].setSVNum(fragment[i].getSVNum());
	       bind[bindingCount].insertSite(i-start);
	       bindingCount++;
	    }
	 }
      }
   }
}

void initialMapping(int clusterNum)
{
   int num = cluster[clusterNum].getFragNum();

   if(bindingCount == 0)
   {
      mapping2geno = new int [num];
      mapping2pheno = new int [num];
      for(int i=0;i<num;i++)
      {
	 mapping2geno[i] = i;
	 mapping2pheno[i] = i;
      }
   }
   else
   {
      int constrainNum = 0, curNum = bindingCount;
      for(int i=0;i<bindingCount;i++) constrainNum = constrainNum + (bind[i].getSize() - 1);
      mapping2geno = new int [num-constrainNum];
      mapping2pheno = new int [num];
      for(int i=0;i<num;i++) mapping2pheno[i] = -1;
      for(int i=0;i<num-constrainNum;i++) mapping2geno[i] = -1;
      for(int i=0;i<bindingCount;i++)
      {
	 for(int j=0;j<bind[i].getSize();j++)
	 {
	    mapping2pheno[bind[i].getSite(j)] = i;
	 }
      }
      for(int i=0;i<num;i++)
      {
	 if(mapping2pheno[i] == -1) 
	 {
	    mapping2pheno[i] = curNum;
	    mapping2geno[curNum] = i;
	    curNum++;
	 }
      }
   }
}

void pipeline()
{
   int i, j, k, haploIndex = 0;
   int *parent = new int [2];
   bool previous = false;
   ofstream fout;
   fout.open("haplotype_validation.log");

   haplotype[0] = new char [snpCount+clusterCount];
   haplotype[1] = new char [snpCount+clusterCount];
   haplotype[0][snpCount+clusterCount-1] = '\0';
   haplotype[1][snpCount+clusterCount-1] = '\0';
   svList = new int* [svCount];
   for(i=0;i<svCount;i++) svList[i] = new int [2];

   for(i=0;i<clusterCount;i++)
   {
      if(cluster[i].getFragNum() == 0)
      {
	 if(previous == true)
	 {
	    haplotype[0][haploIndex] = '|';
	    haplotype[1][haploIndex] = '|';
	    haploIndex++;
	 }
	 for(j=0;j<cluster[i].getSize();j++)
	 {
	    haplotype[0][haploIndex] = '-';
	    haplotype[1][haploIndex] = '-';
	    haploIndex++;
	 }
	 previous = true;
	 continue;
      }
      checkBinding(i);
      initialMapping(i);
      population = new CHROMOSOME [popSize];
      children = new CHROMOSOME [popSize];
      for(j=0;j<popSize;j++)
      {
	 int fragStart = cluster[i].getFragStart();
	 int fragNum = cluster[i].getFragNum();
	 population[j].insertInfo(&fragment[fragStart],bind,bindingCount,fragNum,mapping2pheno,mapping2geno,svList,&svNum);
      }
      thrd = new pthread_t [popSize];
      for(j=0;j<popSize;j++) pthread_create(&thrd[j],NULL,parallel_init,(void*)j);
      for(j=0;j<popSize;j++) pthread_join(thrd[j],NULL);
      for(j=0;j<generation;j++)
      {
	 for(k=0;k<popSize;k+=2)
	 {
	    parentSelection(parent);
	    crossOver(parent,k);
	 }
	 for(k=0;k<popSize;k++) pthread_create(&thrd[k],NULL,parallel,(void*)k);
	 for(k=0;k<popSize;k++) pthread_join(thrd[k],NULL);
	 sortPop();
	 sortChild();
	 survival();
	 sortPop();
      }
      if(previous == true)
      {
	 haplotype[0][haploIndex] = '|';
	 haplotype[1][haploIndex] = '|';
	 haploIndex++;
      }
      int len = fragment[cluster[i].getFragStart()].getSNPNum();
      population[0].updateSVHaplotype();
      for(j=0;j<len;j++)
      {
	 haplotype[0][haploIndex] = population[0].getHaplotype(0)[j];
	 haplotype[1][haploIndex] = population[0].getHaplotype(1)[j];
	 haploIndex++;
      }
      previous = true;
      delete [] bind;
      delete [] mapping2geno;
      delete [] mapping2pheno;
      delete [] population;
      delete [] children;
      bind = NULL;
      population = NULL;
      children = NULL;
      bindingCount = 0;
   }
   fout << haplotype[0] << endl << haplotype[1] << endl;
   fout.close();
   printHaplotype();
   delete [] haplotype[0];
   delete [] haplotype[1];
   delete [] parent;
}

void parentSelection(int *parent)
{
   for(int i=0;i<2;i++)
   {
      int randIndex1 = rand()%popSize;
      int randIndex2 = rand()%popSize;
      if(population[randIndex1].getErrorCorrectionNum() < population[randIndex2].getErrorCorrectionNum()) parent[i] = randIndex1;
      else parent[i] = randIndex2;
   }
}

void crossOver(int *parent, int index)
{
   int size = population[0].getChrSize();

   children[index] = population[parent[0]];
   children[index+1] = population[parent[1]];

   for(int i=0;i<size;i++)
   {
      if(rand()%2 == 1)
      {
	 children[index].setGene(i,population[parent[1]].getGene(i));
	 children[index+1].setGene(i,population[parent[0]].getGene(i));
      }
   }
}

void mutation(CHROMOSOME* pop)
{
   for(int i=0;i<pop->getErrorListSize();i++)
   {
      if(rand()%100 < mutateRate)
      {
	 pop->flipGene(pop->getErrorSite(i));
      }
   }
   int size = pop->getChrSize();
   for(int i=0;i<pop->getNormalListSize();i++)
   {
      if(rand()%size == 0)
      {
	 pop->flipGene(pop->getNormalSite(i));
      }
   }
}

void sortPop()
{
   qsort(population,popSize,sizeof(class CHROMOSOME),compareErrorNum);
}

void sortChild()
{
   qsort(children,popSize,sizeof(class CHROMOSOME),compareErrorNum);
}

void survival()
{
   int popTop = 0, popBot = popSize - 1, childTop = 0;

   while(popTop <= popBot)
   {
      if(population[popTop].getErrorCorrectionNum() <= children[childTop].getErrorCorrectionNum()) popTop++;
      else
      {
	 swapChromosome(population[popBot],children[childTop]);
	 popBot--;
	 childTop++;
      }
   }
}

int compareErrorNum(const void *s1, const void *s2)
{
   return (*(class CHROMOSOME*)s1).getErrorCorrectionNum() - (*(class CHROMOSOME*)s2).getErrorCorrectionNum();
}

void* parallel(void* data)
{
   long index = (long)data;
   children[index].evaluate();
   mutation(&children[index]);
   children[index].evaluate();
   return (void*)0;
}

void* parallel_init(void *data)
{
   long index = (long)data;
   population[index].initial();
   if(!population[index].checkChr())
   {
      printf("ERROR! There are some noise in chr.\n");
      printf("%s\n",population[index].getChr());
      exit(0);
   }
   population[index].evaluate();
   return (void*)0;
}

void printHaplotype()
{
   int snp_index = 0, sv_index = 0, previous = 0, size, pos, type;
   FILE *out, *snp, *sv;
   int *SNP = new int [snpCount];
   char ref[1024], buf[1024], *tmp;

   out = fopen("haplotype.log","w");
   snp = fopen("snp_site.log","r");

   for(int i=0;i<snpCount;i++)
   {
      fgets(buf,1024,snp);
      tmp = strtok(buf,"\t");
      strcpy(ref,tmp);
      tmp = strtok(NULL,"\t");
      SNP[i] = atoi(tmp);
   }

   if(svNum != 0) 
   {
      sv = fopen("detected_sv.log","r");
      for(int i=0;i<=svList[sv_index][0];i++) fgets(buf,1024,sv);
      tmp = strtok(buf,"\t");
      if(strcmp(tmp,"Deletion") == 0) type = 1;
      else if(strcmp(tmp,"Insertion") == 0) type = 2;
      else if(strcmp(tmp,"Inversion") == 0) type = 3;
      size = atoi(strtok(NULL,"\t"));
      pos = atoi(strtok(NULL,"\t"));
      previous = svList[sv_index][0];
   }

   fprintf(out,"#ref_name\ttype\tpos\tsize\thap0\thap1\n");

   for(int i=0;i<snpCount+clusterCount-1;i++)
   {
      if(haplotype[0][i] == '|') 
      {
	 fprintf(out,"-\n");
	 continue;
      }
      if(SNP[snp_index] > pos && sv_index < svNum)
      {
	 if(type == 1) fprintf(out,"%s\t%s\t",ref,"Del");
	 else if(type == 2) fprintf(out,"%s\t%s\t",ref,"Ins");
	 else if(type == 3) fprintf(out,"%s\t%s\t",ref,"Rev");
	 fprintf(out,"%d\t%d\t",pos,size);
	 if(svList[sv_index][1] == 0) fprintf(out,"V\t-\n");
	 else fprintf(out,"-\tV\n");
	 sv_index++;
	 if(sv_index < svNum)
	 {
	    for(int i=previous;i<svList[sv_index][0];i++) fgets(buf,1024,sv);
	    tmp = strtok(buf,"\t");
	    if(strcmp(tmp,"Deletion") == 0) type = 1;
	    else if(strcmp(tmp,"Insertion") == 0) type = 2;
	    else if(strcmp(tmp,"Inversion") == 0) type = 3;
	    size = atoi(strtok(NULL,"\t"));
	    pos = atoi(strtok(NULL,"\t"));
	    previous = svList[sv_index][0];
	 }
      }
      fprintf(out,"%s\t%s\t%d\t-\t%c\t%c\n",ref,"SNP",SNP[snp_index],haplotype[0][i],haplotype[1][i]);
      snp_index++;
   }

   fclose(out);
   fclose(snp);
   if(svNum != 0) fclose(sv);
}
