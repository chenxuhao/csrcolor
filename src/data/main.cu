// Copyright 2016, National University of Defense Technology
// Authors: Xuhao Chen <cxh@illinois.edu> and Pingfan Li <lipingfan@163.com>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include "lonestargpu.h"
#include "variants.h"
using namespace std;

#ifndef	ITERATIONS
#define	ITERATIONS 1
#endif
#ifndef	BLKSIZE
#define	BLKSIZE 128
#endif

// transfer R-MAT generated gr graph to CSR format
void gr2csr(char *gr, int &m, int &nnz, int *&csrRowPtr, int *&csrColInd) {
	printf("Reading RMAT (.gr) input file %s\n", gr);
	std::ifstream cfile;
	cfile.open(gr);
	std::string str;
	getline(cfile, str);
	char c;
	sscanf(str.c_str(), "%c", &c);
	while (c == 'c') {
		getline(cfile, str);
		sscanf(str.c_str(), "%c", &c);
	}
	char sp[3];
	sscanf(str.c_str(), "%c %s %d %d", &c, sp, &m, &nnz);
	printf("num_vertices %d num_edges %d\n", m, nnz);
	//printf("%c %s %d %d\n", c, sp, m, nnz);
	vector<set<int> > svector;
	set<int> s;
	for (int i = 0; i < m; i++)
		svector.push_back(s);
	int dst, src;
	for (int i = 0; i < nnz; i++) {
		getline(cfile, str);
		sscanf(str.c_str(), "%c %d %d", &c, &src, &dst);

		if (c != 'a')
			printf("line %d\n", __LINE__);
		dst--;
		src--;
		svector[src].insert(dst);
		svector[dst].insert(src);
	}
	csrRowPtr = (int *)malloc((m + 1) * sizeof(int));
	int count = 0;
	for (int i = 0; i < m; i++) {
		csrRowPtr[i] = count;
		count += svector[i].size();
	}
	csrRowPtr[m] = count;
	if (count != nnz) {
		printf("This graph is not symmetric\n");
		nnz = count;
	}
	double avgdeg;
	double variance = 0.0;
	int maxdeg = 0;
	int mindeg = m;
	avgdeg = (double)nnz / m;
	for (int i = 0; i < m; i++) {
		int deg_i = csrRowPtr[i + 1] - csrRowPtr[i];
		if (deg_i > maxdeg)
			maxdeg = deg_i;
		if (deg_i < mindeg)
			mindeg = deg_i;
		variance += (deg_i - avgdeg) * (deg_i - avgdeg) / m;
	}
	printf("mindeg %d maxdeg %d avgdeg %.2f variance %.2f\n", mindeg, maxdeg, avgdeg, variance);
	csrColInd = (int *)malloc(count * sizeof(int));
	set<int>::iterator site;
	for (int i = 0, index = 0; i < m; i++) {
		site = svector[i].begin();
		while (site != svector[i].end()) {
			csrColInd[index++] = *site;
			site++;
		}
	}
}

// transfer *.graph file to CSR format
void graph2csr(char *graph, int &m, int &nnz, int *&csrRowPtr, int *&csrColInd) {
	printf("Reading .graph input file %s\n", graph);
	std::ifstream cfile;
	cfile.open(graph);
	std::string str;
	getline(cfile, str);
	sscanf(str.c_str(), "%d %d", &m, &nnz);
	printf("num_vertices %d num_edges %d\n", m, nnz);
	vector<set<int> > svector;
	set<int> s;
	for (int i = 0; i < m; i++)
		svector.push_back(s);
	int dst;
	for (int i = 0; i < m; i++) {
		getline(cfile, str);
		istringstream istr;
		istr.str(str);
		while(istr>>dst) {
			dst --;
			svector[i].insert(dst);
			svector[dst].insert(i);
		}
		istr.clear();
	}
    cfile.close();
	csrRowPtr = (int *)malloc((m + 1) * sizeof(int));
	int count = 0;
	for (int i = 0; i < m; i++) {
		csrRowPtr[i] = count;
		count += svector[i].size();
	}
	csrRowPtr[m] = count;
	if (count != nnz) {
		printf("This graph is not symmetric\n");
		nnz = count;
	}
	double avgdeg;
	double variance = 0.0;
	int maxdeg = 0;
	int mindeg = m;
	avgdeg = (double)nnz / m;
	for (int i = 0; i < m; i++) {
		int deg_i = csrRowPtr[i + 1] - csrRowPtr[i];
		if (deg_i > maxdeg)
			maxdeg = deg_i;
		if (deg_i < mindeg)
			mindeg = deg_i;
		variance += (deg_i - avgdeg) * (deg_i - avgdeg) / m;
	}
	printf("mindeg %d maxdeg %d avgdeg %.2f variance %.2f\n", mindeg, maxdeg, avgdeg, variance);
	csrColInd = (int *)malloc(count * sizeof(int));
	set<int>::iterator site;
	for (int i = 0, index = 0; i < m; i++) {
		site = svector[i].begin();
		while (site != svector[i].end()) {
			csrColInd[index++] = *site;
			site++;
		}
	}
}

// transfer mtx graph to CSR format
void mtx2csr(char *mtx, int &m, int &nnz, int *&csrRowPtr, int *&csrColInd) {
	printf("Reading (.mtx) input file %s\n", mtx);
	std::ifstream cfile;
	cfile.open(mtx);
	std::string str;
	getline(cfile, str);
	char c;
	sscanf(str.c_str(), "%c", &c);
	while (c == '%') {
		getline(cfile, str);
		sscanf(str.c_str(), "%c", &c);
	}
	int n;
	sscanf(str.c_str(), "%d %d %d", &m, &n, &nnz);
	if (m != n) {
		printf("error!\n");
		exit(0);
	}
	printf("num_vertices %d num_edges %d\n", m, nnz);
	vector<set<int> > svector;
	set<int> s;
	for (int i = 0; i < m; i++)
		svector.push_back(s);
	int dst, src;
	for (int i = 0; i < nnz; i++) {
		getline(cfile, str);
		sscanf(str.c_str(), "%d %d", &dst, &src);
		dst--;
		src--;
		svector[src].insert(dst);
		svector[dst].insert(src);
	}
	cfile.close();
	csrRowPtr = (int *)malloc((m + 1) * sizeof(int));
	int count = 0;
	for (int i = 0; i < m; i++) {
		csrRowPtr[i] = count;
		count += svector[i].size();
	}
	csrRowPtr[m] = count;
	if (count != nnz) {
		printf("This graph is not symmetric\n");
		nnz = count;
	}
	double avgdeg;
	double variance = 0.0;
	int maxdeg = 0;
	int mindeg = m;
	avgdeg = (double)nnz / m;
	for (int i = 0; i < m; i++) {
		int deg_i = csrRowPtr[i + 1] - csrRowPtr[i];
		if (deg_i > maxdeg)
			maxdeg = deg_i;
		if (deg_i < mindeg)
			mindeg = deg_i;
		variance += (deg_i - avgdeg) * (deg_i - avgdeg) / m;
	}
	printf("mindeg %d maxdeg %d avgdeg %.2f variance %.2f\n", mindeg, maxdeg, avgdeg, variance);
	csrColInd = (int *)malloc(count * sizeof(int));
	set<int>::iterator site;
	for (int i = 0, index = 0; i < m; i++) {
		site = svector[i].begin();
		while (site != svector[i].end()) {
			csrColInd[index++] = *site;
			site++;
		}
	}
}

// store colour of all vertex
void write_solution(char *fname, int *coloring, int n) {
	int i;
	FILE *fp;
	fp = fopen(fname, "w");
	for (i = 0; i < n; i++) {
		//fprintf(fp, "%d:%d\n", i, coloring[i]);
		fprintf(fp, "%d\n", coloring[i]);
	}
	fclose(fp);
}

// check if correctly coloured
void verify(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int *correct) {
	int i, offset, neighbor_j;
	for (i = 0; i < m; i++) {
		for (offset = csrRowPtr[i]; offset < csrRowPtr[i + 1]; offset++) {
			neighbor_j = csrColInd[offset];
			if (coloring[i] == coloring[neighbor_j] && neighbor_j != i) {
				*correct = 0;
				//printf("coloring[%d] = coloring[%d] = %d\n", i, neighbor_j, coloring[i]);
				break;
			}
		}	
	}
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: %s <graph> <num_SMs>\n", argv[0]);
		exit(1);
	}
	int m, nnz, *csrRowPtr = NULL, *csrColInd = NULL;
	// read graph
	if (strstr(argv[1], ".mtx"))
		mtx2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[1], ".graph"))
		graph2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[1], ".gr"))
		gr2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else { printf("Unrecognizable input file format\n"); exit(0); }
	int *coloring = (int *)calloc(m, sizeof(int));
	int correct = 1;
	int num_SMs;
	if (argc > 2) {
		num_SMs = atoi(argv[2]);
		printf("block_size=%d, num_SMs=%d\n", BLKSIZE, num_SMs);
	}
#if VARIANT==DATA_LDB
	color_ldb(m, nnz, csrRowPtr, csrColInd, coloring, num_SMs);
#else
	color(m, nnz, csrRowPtr, csrColInd, coloring, num_SMs);
#endif
	write_solution("color.txt", coloring, m);
	verify(m, nnz, csrRowPtr, csrColInd, coloring, &correct);
	if (correct)
		printf("correct.\n");
	else
		printf("incorrect.\n");
	return 0;
}
