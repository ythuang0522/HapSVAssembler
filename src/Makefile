OPT = -O2
CC = g++

all: snpMatrixCreator genomeCreator svDetector hapAssembler

hapAssembler: cluster.o fragment.o bindingSite.o chromosome.o hapAssembler.cpp
	$(CC) -o ../bin/hapAssembler cluster.o fragment.o bindingSite.o chromosome.o hapAssembler.cpp -pthread $(OPT)

svDetector: svDetector.cpp
	$(CC) -o ../bin/svDetector svDetector.cpp $(OPT)

snpMatrixCreator: snp.o sv.o snpMatrixCreator.cpp
	$(CC) -o ../bin/snpMatrixCreator snp.o sv.o snpMatrixCreator.cpp $(OPT)

genomeCreator: genomeCreator.cpp
	$(CC) -o ../bin/genomeCreator genomeCreator.cpp $(OPT)

cluster.o: cluster.cpp
	$(CC) -c cluster.cpp $(OPT)

fragment.o: fragment.cpp
	$(CC) -c fragment.cpp $(OPT)

bindingSite.o: bindingSite.cpp
	$(CC) -c bindingSite.cpp $(OPT)

chromosome.o: chromosome.cpp
	$(CC) -c chromosome.cpp $(OPT)

snp.o: snp.cpp
	$(CC) -c snp.cpp $(OPT)

sv.o: sv.cpp
	$(CC) -c sv.cpp $(OPT)

clean_all:
	rm -f *.o ../bin/svDetector ../bin/snpMatrixCreator ../bin/hapAssembler \
      ../haplotype_result/* ../bwa_result/* ../denovo_result/* 

clean_haplotype:
	rm -f *.o ../bin/svDetector ../bin/snpMatrixCreator ../bin/hapAssembler \
      ../haplotype_result/*
