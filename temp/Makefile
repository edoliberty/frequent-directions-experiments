main: common.c sparseVector.c sparseMatrix.c sparseSketcher.c usage.c
	gcc -c -g -o common.o common.c
	gcc -c -g -o sparseVector.o sparseVector.c
	gcc -c -g -o sparseMatrix.o sparseMatrix.c	
	gcc -c -g -o sparseSketcher.o sparseSketcher.c		
	gcc -c -g -o usage.o usage.c	
	gfortran common.o sparseVector.o sparseMatrix.o sparseSketcher.o usage.o /usr/lib/liblapacke.a /usr/lib/liblapack.a /usr/lib/libblas.a -o usage.exe
