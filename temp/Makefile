
main: 
	gcc -c -O3 -g -std=c99 -o sparseVector.o sparseVector.c
	gcc -c -O3 -g -std=c99 -o sparseMatrix.o sparseMatrix.c
	gcc -c -O3 -g -std=c99 -o sparseSketcher.o sparseSketcher.c
	gcc -c -O3 -g -std=c99 -o frequentDirections.o frequentDirections.c	
	gcc -c -O3 -g -std=c99 -o testAll.o testAll.c
	gcc -c -O3 -g -std=c99 -o common.o common.c
	gfortran common.o frequentDirections.o sparseVector.o sparseMatrix.o sparseSketcher.o testAll.o /usr/lib/liblapacke.a /usr/lib/liblapack.a /usr/lib/libblas.a -o testAll.exe


