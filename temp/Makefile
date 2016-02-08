main: common.c sparseVector.c sparseMatrix.c sparseSketcher.c usage.c
	gcc -c -g -o common.o common.c
	gcc -c -g -o sparseVector.o sparseVector.c
	gcc -c -g -o sparseMatrix.o sparseMatrix.c	
	gcc -c -g -o sparseSketcher.o sparseSketcher.c		
	gcc -c -g -o usage.o usage.c	
	gfortran common.o sparseVector.o sparseMatrix.o sparseSketcher.o usage.o /usr/lib/liblapacke.a /usr/lib/liblapack.a /usr/lib/libblas.a -o usage.exe


randomSum: randomSum.c common.c usage.c sparseVector.c
	gcc -c -o common.o common.c
	gcc -c -o sparseVector.o sparseVector.c
	gcc -c -o randomSum.o randomSum.c
	gcc -c -o randomProjection.o randomProjection.c
	gcc -c -g -o sparseMatrix.o sparseMatrix.c	
	gcc -c -o sparseSketcher.o sparseSketcher.c
	gcc -c -o usage.o usage.c
	gfortran common.o sparseVector.o randomSum.o randomProjection.o sparseMatrix.o sparseSketcher.o usage.o /usr/lib/liblapacke.a /usr/lib/liblapack.a /usr/lib/libblas.a -o usage.exe

itemsampler: sparseVector.c singleItemSampler.c
	gcc -c -o sparseVector.o sparseVector.c
	gcc -c -o singleItemSampler.o singleItemSampler.c
	gcc sparseVector.o singleItemSampler.o -lm -o merge.o


rowsampler: sparseVector.c singleItemSampler.c rowSampler.c
	gcc -c -o sparseVector.o sparseVector.c
	gcc -c -o singleItemSampler.o singleItemSampler.c
	gcc -c -o rowSampler.o rowSampler.c	
	gcc -c -o usage.o usage.c
	gcc -c -o common.o common.c
	gcc common.o sparseVector.o singleItemSampler.o rowSampler.o usage.o -lm -o merge.o


all: 
	gcc -c -g -o common.o common.c
	gcc -c -o singleItemSampler.o singleItemSampler.c
	gcc -c -o rowSampler.o rowSampler.c	
	gcc -c -o randomSum.o randomSum.c
	gcc -c -o randomProjection.o randomProjection.c
	gcc -c -g -o sparseMatrix.o sparseMatrix.c	
	gcc -c -g -o sparseVector.o sparseVector.c
	gcc -c -g -o frequentDirections.o frequentDirections.c
	gcc -c -g -o sparseSketcher.o sparseSketcher.c		
	gcc -c -g -o testAll.o testAll.c	
	gfortran common.o frequentDirections.o rowSampler.o randomSum.o randomProjection.o singleItemSampler.o sparseVector.o sparseMatrix.o sparseSketcher.o testAll.o /usr/lib/liblapacke.a /usr/lib/liblapack.a /usr/lib/libblas.a -o testAll.exe
