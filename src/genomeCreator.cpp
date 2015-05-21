#include <unistd.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
using namespace std;

struct SV
{
	int type;
	int size;
	int leftPos;
	int rightPos;
};

struct SV *svList;
char *haplo1, *haplo2;
char *genome1_with_N, *genome1_without_N, *genome2_without_N, *genome2_with_N, *refOrigin;
char *string_n, *string_insert;
bool *snp;
int *insertList;
char *tmpList;
int svBase = 300;

char fastaName[64];
int genomeStart = 0;
int genomeSize = 10000;
int readLength = 75;
int coverage = 20;
int difference = 1;
long readCount = 0;
long readCount1 = 0, readCount2 = 0;
int snpRate = 100;
int svCount = 0;
int deleteCount = 0;
int insertCount = 0;
int reverseCount = 0;
int svFlag = 0;
int window = 0;
int totInsertSize = 0;

static int usage();
void memoryAllocate();
void decideCount();
void countReadNum();
void createRef(char *name);
void selectPos();
void permutation();
void createDeletion(int index);
void createInsertion(int index);
void createReversion(int index);
void snpMask();
void createSNP();
void outputHaplotype();
void outputGenome();
void memoryRecycle();

int main(int argc, char *argv[])
{
   char c;
   strcpy(fastaName,"");

   while((c = getopt(argc,argv,"f:S:g:c:l:s:D:dirh")) >= 0)
   {
		switch(c)
		{
			case 'f': strcpy(fastaName,optarg); break;
			case 'S': genomeStart = atoi(optarg); break;
			case 'g': genomeSize = atoi(optarg); break;
			case 'c': coverage = atoi(optarg); break;
			case 'l': readLength = atoi(optarg); break;
			case 's': snpRate = atoi(optarg); break;
			case 'd': svFlag = svFlag | 4; break;
			case 'i': svFlag = svFlag | 2; break;
			case 'r': svFlag = svFlag | 1; break;
			case 'D': difference = atoi(optarg); break;
			case 'h': usage(); return 0;
			default : usage(); return 1;
		}
	}
	if(optind < 1) return usage();
	if(snpRate < 100 || snpRate > 10000) 
	{
	   	printf("SNP rate s should be in [1,10000], input argument s = %d\n\n",snpRate);
		exit(1);
	}
	if(difference < 1 || difference > 20)
	{
	   	printf("difference rate D should be in [1,5], input argument D = %d\n\n",difference);
	   	exit(1);
	}

	double startT, endT;
	double overallTime = 0, sectionTime = 0;

	startT = clock();   
	srand(time(NULL));
	readCount = ((float)coverage*(float)genomeSize)/(2*readLength);
	decideCount();
	memoryAllocate();
	createRef(fastaName);
	printf("[createRef] genome size: %d\n",genomeSize);
	
	if(svCount > 0) 
	{
	   	selectPos();
	   	printf("[createSV] summary of SVs: ");
	  	if(deleteCount > 0) printf("%d deletions ",deleteCount);
	   	if(insertCount > 0) printf("%d insertions ",insertCount);
	   	if(reverseCount > 0) printf("%d inversions ",reverseCount);
		permutation();
		printf("\n[createSV] SVs have been generated.\n");
	}

	snpMask();
	createSNP();
	printf("[createSNP] SNPs have been generated. There are %ld SNP sites.\n", strlen(haplo1));
	outputHaplotype();
	countReadNum();
	outputGenome();
	memoryRecycle();
	endT = clock();
	printf("[countReadNum] number of reads from chromatid1: %ld\n",readCount1);
	printf("[countReadNum] number of reads from chromatid2: %ld\n",readCount2);
	printf("[genomeCreator] time elapses: %.2lf Seconds\n",(endT-startT)/CLOCKS_PER_SEC);
	
	return 0;	
}

void memoryRecycle()
{
	delete [] haplo1;
	delete [] haplo2;
	delete [] refOrigin;
	delete [] genome1_with_N;
	delete [] genome1_without_N;
	delete [] genome2_with_N;
	delete [] genome2_without_N;
	delete [] svList;
	delete [] insertList;
	delete [] string_n;
	delete [] snp;
}

void decideCount()
{
//   	svCount = (float)difference/100*genomeSize/300;
   	svCount = (float)genomeSize/100/300;
	(svCount==0)?(svCount=1):(svCount=svCount);

	switch(svFlag)
	{
	   	case 0: svCount = 0;
		case 1: reverseCount = svCount; break;
		case 2: insertCount = svCount; break;
		case 4: deleteCount = svCount; break;
		case 3: reverseCount = svCount/2; insertCount = svCount - reverseCount; break;
		case 5: reverseCount = svCount/2; deleteCount = svCount - reverseCount; break;
		case 6: deleteCount = svCount/2; insertCount = svCount - deleteCount; break;
		case 7: deleteCount = svCount/3; reverseCount = deleteCount; insertCount = svCount - deleteCount - reverseCount; break;
		default: break;
	}
}

void memoryAllocate()
{
	haplo1 = new char [genomeSize];
	haplo2 = new char [genomeSize];
	genome1_with_N= new char [2*genomeSize];
	genome1_without_N = new char [2*genomeSize];
	genome2_with_N = new char [2*genomeSize];
	genome2_without_N = new char [2*genomeSize];
	snp = new bool [2*genomeSize];
	tmpList = new char [2*genomeSize];
}

void outputGenome()
{
	ofstream fout1, fout2, fout3;
	int index = 0, i = 0;

	fout1.open("reference.fasta");
	fout2.open("reference_chromatid2.fasta");
	fout3.open("aligned_reference.fasta");

	while(genome1_with_N[i] != '\0')
	{
		if(genome1_with_N[i] != 'N')
		{
			genome1_without_N[index] = genome1_with_N[i];
			index++;
		}
		i++;
	}
	genome1_without_N[index] = '\0';
	index = 0;
	i = 0;
	while(genome2_with_N[i] != '\0')
	{
		if(genome2_with_N[i] != 'N')
		{
			genome2_without_N[index] = genome2_with_N[i];
			index++;
		}
		i++;
	}
	genome2_without_N[index] = '\0';

	fout1 << ">reference\n" << genome1_without_N << endl;
	fout2 << ">reference_chromatid2\n" << genome2_without_N << endl;
	fout3 << ">reference_1" << endl << genome1_with_N << endl;
	fout3 << ">reference_2" << endl << genome2_with_N << endl;
	fout1.close();
	fout2.close();
	fout3.close();
}

void countReadNum()
{
	for(int i=0;i<readCount;i++)
	{
		if(rand()%2 == 0) readCount1++;
		else readCount2++;
	}
}

void createRef(char *name)
{
	int allele;

	if(strcmp(fastaName,"") == 0)
	{
	   	for(int i=0;i<genomeSize;i++)
		{
			allele = rand()%4;
			if(allele == 0) genome1_with_N[i] = 'A';
			else if(allele == 1) genome1_with_N[i] = 'C';
			else if(allele == 2) genome1_with_N[i] = 'T';
			else genome1_with_N[i] = 'G';
		}
		genome1_with_N[genomeSize] = '\0';
	}
	else
	{
	   	ifstream fin;
		int len = 0;
		char buf[1024], c;

		fin.open(name);
		fin.getline(buf,1024);
		if(buf[0] != '>') 
		{
		   	cout << "[ERROR] FASTA file format error!\n\n";
			exit(0);
		}
		while((c = fin.get()) != EOF)
		{
		   	if(c != '\n') len++;
		}
		refOrigin = new char [len+1];
		int index = 0;
		fin.close();
		fin.open(name);
		fin.getline(buf,1024);
		while((c = fin.get()) != EOF)
		{
		   	if(c != '\n')
			{
			   	if(toupper(c) == 'N')
				{
				   	allele = rand()%4;
					if(allele == 0) refOrigin[index] = 'A';
					else if(allele == 1) refOrigin[index] = 'C';
					else if(allele == 2) refOrigin[index] = 'T';
					else refOrigin[index] = 'G';
				}
				else refOrigin[index] = toupper(c);
				index++;
			}
		}
		fin.close();
		refOrigin[index] = '\0';
		if((genomeStart+genomeSize) > len)
		{
		   	printf("[ERROR] The bound is out of range(%d,%d)!\n\n",0,len-1);
			exit(0);
		}
		for(int i=0;i<genomeSize;i++)
		{
		   	c = refOrigin[i+genomeStart];
			genome1_with_N[i] = c;
		}
	}
	strcpy(genome2_with_N,genome1_with_N);
}

void selectPos()
{
	window = genomeSize / svCount;
	int i = 0;
	int max = 0;

	svList = new class SV [svCount];

	for(int j=0;j<deleteCount;j++)
	{
		svList[i].type = 0;
//		svList[i].size = rand()%(window/50 - 50) + 50;
		svList[i].size = rand()%400 + (300*difference-200);
		svList[i].leftPos = rand()%(window-svList[i].size);
		svList[i].rightPos = svList[i].leftPos + svList[i].size - 1;
		i++;
	}
	for(int j=0;j<insertCount;j++)
	{
		svList[i].type = 1;
//		svList[i].size = rand()%(window/50 - 50) + 50;
		svList[i].size = rand()%400 + (300*difference-200);
		if(svList[i].size > max) max = svList[i].size;
		svList[i].leftPos = rand()%window;
		svList[i].rightPos = svList[i].leftPos + 1;
		i++;
	}
	for(int j=0;j<reverseCount;j++)
	{
		svList[i].type = 2;
//		svList[i].size = rand()%(window/20 - 150) + 150;
		svList[i].size = rand()%400 + (300*difference-200);
		svList[i].leftPos = rand()%(window-svList[i].size);
		svList[i].rightPos = svList[i].leftPos + svList[i].size - 1;
		i++;
	}
	string_n = new char [max+1];
	string_insert = new char [max+1];
	for(int j=0;j<max;j++) string_n[j] = 'N';
	string_n[max] = '\0';
}	

void permutation()
{
	ofstream fout;
	struct SV tmp;
	int e1, e2, index = 0;
	insertList = new int [insertCount];

	fout.open("structural_variation.log");

	if(deleteCount!=svCount && insertCount!=svCount && reverseCount!=svCount)
	{
	   	for(int i=0;i<svCount/2;i++)
		{
			e1 = rand()%svCount;
			e2 = rand()%svCount;
			tmp = svList[e1];
			svList[e1] = svList[e2];
			svList[e2] = tmp;
		}
	}
	for(int i=0;i<svCount;i++)
	{
	   	switch(svList[i].type)
		{
			case 0: fout << "Deletion\t"; createDeletion(i); break;
			case 1: fout << "Insertion\t"; insertList[index] = i; index++; break;
			case 2: fout << "Inversion\t"; createReversion(i); break;
		}
		fout << svList[i].size << "\t" << i*window+svList[i].leftPos << "\t" << i*window+svList[i].rightPos << endl;
	}
	
	fout.close();
	for(int i=0;i<index;i++) createInsertion(insertList[i]);
}

void createDeletion(int index)
{
	for(int i=0;i<svList[index].size;i++) genome2_with_N[index*window+svList[index].leftPos+i] = 'N';
}

void createInsertion(int index)
{
	int allele, i;
	int offset = index*window + totInsertSize;
	int rightBP = offset+svList[index].rightPos;

	strcpy(tmpList,&genome1_with_N[rightBP]);
	strncpy(&genome1_with_N[rightBP],string_n,svList[index].size);
	strcpy(&genome1_with_N[rightBP+svList[index].size],tmpList);

	strcpy(tmpList,&genome2_with_N[rightBP]);
	for(i=0;i<svList[index].size;i++)
	{
	   	allele = rand()%4;
		if(allele == 0) string_insert[i] = 'A';
		else if(allele == 1) string_insert[i] = 'C';
		else if(allele == 2) string_insert[i] = 'T';
		else string_insert[i] = 'G';
	}
	string_insert[i] = '\0';
	strcpy(&genome2_with_N[rightBP],string_insert);
	strcpy(&genome2_with_N[rightBP+svList[index].size],tmpList);
	totInsertSize = totInsertSize + svList[index].size;
}

void createReversion(int index)
{
	char allele;

	for(int i=0;i<svList[index].size;i++)
	{
		switch(genome1_with_N[index*window+svList[index].leftPos+svList[index].size-i])
		{
			case 'A': allele = 'T'; break;
			case 'C': allele = 'G'; break;
			case 'T': allele = 'A'; break;
			case 'G': allele = 'C'; break;
		}
		genome2_with_N[index*window+svList[index].leftPos+i] = allele;
	}
}

void snpMask()
{
   	int offset = 0;

	for(int i=0;i<2*genomeSize;i++) snp[i] = true;
   	for(int i=0;i<svCount;i++)
	{
	   	if(svList[i].type == 0 || svList[i].type == 2) //deletion or inversion
		{
			for(int j=svList[i].leftPos;j<=svList[i].rightPos;j++) snp[offset+j] = false;
		}
		else //insertion
		{
			offset = offset + svList[i].size;
		}
		offset = offset + window;
	}
}

void createSNP()
{
	ofstream fout;
	int index = 0, ref_index = 0;
	int length = strlen(genome1_with_N);
	fout.open("snp_site.log");
	for(int i=0;i<length;i++)
	{
		if(genome1_with_N[i] == 'N') continue;
		if(genome1_with_N[i] != 'N') ref_index++;
	   	if(snp[i] == false) continue;
		int tmp = rand()%snpRate;
		if(tmp == 0)
		{
		   	tmp = rand()%4;
			if(tmp == 0) genome2_with_N[i] = 'A';
			else if (tmp == 1) genome2_with_N[i] = 'C';
			else if (tmp == 2) genome2_with_N[i] = 'T';
			else if (tmp == 3) genome2_with_N[i] = 'G';
			if(genome1_with_N[i]!=genome2_with_N[i]) 
			{
				haplo1[index] = genome1_with_N[i];
				haplo2[index] = genome2_with_N[i];
				index++;
				fout << "reference\t" << ref_index-1 << "\t" << genome1_with_N[i] << "\t" << genome2_with_N[i] << endl;
			}
		}
	}
	haplo1[index] = '\0';
	haplo2[index] = '\0';

	fout.close();
}

void outputHaplotype()
{
	ofstream fout;

	fout.open("haplotype");
	fout << ">haplotype_1\t" << strlen(haplo1) << endl;
	fout << haplo1 << endl;
	fout << ">haplotype_2\t" << strlen(haplo2) << endl;
	fout << haplo2 << endl;
	fout.close();
}

static int usage()
{
	cout << "\nVersion 1.1: released on July 7th, 2011" << endl;
	cout << "\nUsage: genomeCreator [OPTION]" << endl << endl;
	cout << "  [OPTION]" << endl;
	cout << "  -f <FILE> filename of reference genome in FASTA format" << endl;
	cout << "  -S <int>  started-cutting position on reference [Default 0] (required only when -f <FILE> is set)" << endl;
	cout << "  -g <int>  genome size [Default 10000]" << endl;
	cout << "  -c <int>  average coverage of paired-end reads [Default 20]" << endl;
	cout << "  -l <int>  number of known nucleotides on each ends [Default 75]" << endl;
	cout << "  -s <int>  SNP rate (1/s, s=[100,10000]) in diploid genome [Default 100]" << endl;
	cout << "  -D <int>  difference rate (D%, s=[1,10]) between two homologous chromosomes [Default 1]" << endl;
	cout << "  -d        create deletion polymorphism [Defalut NO]" << endl;
	cout << "  -i        create insertion polymorphism [Defalut NO]" << endl;
	cout << "  -r        create reversion(inversion) polymorphism [Defalut NO]" << endl;
	cout << "  -h        help" << endl << endl;
	return 1;
}
