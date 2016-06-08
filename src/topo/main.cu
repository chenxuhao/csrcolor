#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
using namespace std;
#ifndef	ITERATIONS
#define	ITERATIONS 10
#endif
#ifndef	BLKSIZE
#define	BLKSIZE	32
#endif
#include "kernel.h"
#include "graph_io.h"

int main(int argc, char *argv[]) {
	if (argc != 3) {
		printf("Usage: %s <BLKSIZE> <graph>\n", argv[0]);
		exit(1);
	}
	int blksize = atoi(argv[1]);
	int m, nnz, *csrRowPtr = NULL, *csrColInd = NULL;
	if (strstr(argv[2], ".mtx"))
		mtx2csr(argv[2], m, nnz, csrRowPtr, csrColInd);
	if (strstr(argv[2], ".gr"))
		gr2csr(argv[2], m, nnz, csrRowPtr, csrColInd);
	int *coloring = (int *)calloc(m, sizeof(int));
	color(m, nnz, csrRowPtr, csrColInd, coloring, blksize);
	write_solution("color.txt", coloring, m);
	int correct = 1;
	verify(m, nnz, csrRowPtr, csrColInd, coloring, &correct);
	if (correct)
		printf("correct.\n");
	else
		printf("incorrect.\n");
	return 0;
}
