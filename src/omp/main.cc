// Copyright 2016, National University of Defense Technology
// Authors: Xuhao Chen <cxh@illinois.edu> and Pingfan Li <lipingfan@163.com>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <inttypes.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <vector>
#include <set>
#include "common.h"
#include "worklist.h"
typedef unsigned foru;
#include "graph.h"
int num_omp_threads;
using namespace std;
//#include "kernel1.h"
#include "kernel2.h"

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
	nnz = count;
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
		printf("The graph is not symmetric\n");
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

void verify(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int *correct) {
	int i, offset, neighbor_j;
	for (i = 0; i < m; i++) {
		for (offset = csrRowPtr[i]; offset < csrRowPtr[i + 1]; offset++) {
			neighbor_j = csrColInd[offset];
			if (coloring[i] == coloring[neighbor_j] && neighbor_j != i) {
				*correct = 0;
				printf("coloring[%d] = coloring[%d] = %d\n", i, neighbor_j, coloring[i]);
				break;
			}
		}
	}
}

void write_solution(char *fname, int nnodes, int *coloring) {
	int i;
	FILE *fp = fopen(fname, "w");
	for (i = 0; i < nnodes; i++) {
		fprintf(fp, "%d\n", coloring[i]);
	}
	fclose(fp);
}

void mtx2edges(char *mtx, char *edges) {
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
	int m, n, nnz;
	sscanf(str.c_str(), "%d %d %d", &m, &n, &nnz);
	if (m != n) {
		printf("error!\n");
		exit(0);
	}
	vector<set<int> > svector;
	set<int> s;
	for (int i = 0; i < m; i++)
		svector.push_back(s);

	FILE *fp = fopen(edges, "w");

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
	int count = 0;
	for (int i = 0; i < m; i++) {
		count += svector[i].size();
	}
	fprintf(fp, "%d %d\n", m, count);
	set<int>::iterator site;
	for (int i = 0; i < m; i++) {
		site = svector[i].begin();
		while (site != svector[i].end()) { 
			fprintf(fp, "%d %d\n", i, *site);
			site++;
		}
	}
	fclose(fp);
}

void verify(Graph &graph, int *coloring, int *correct) {
	int nnodes = graph.nnodes;
	int i, j, neighbors, neighbor_j;	
	for (i = 0; i < nnodes; i++) {
		neighbors = graph.noutgoing[i];
		for (j = 0; j < neighbors; j++) {
			neighbor_j = graph.edgessrcdst[graph.psrc[i] + j];
			if (coloring[i] == coloring[neighbor_j] && neighbor_j != i) {
				*correct = 0;
				printf("colors[%d] = colors[%d] = %d\n", i, neighbor_j, coloring[i]);
				break;
			}
		}	
	}
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		printf("Usage: %s <nThreads> <graph>\n", argv[0]);
		exit(1);
	}
	num_omp_threads = atoi(argv[1]);
#ifdef ENABLE_OPENMP
	omp_set_num_threads(num_omp_threads);
	printf("OpenMP graph coloring by Xuhao Chen, num_omp_threads=%d\n", num_omp_threads);
#endif
	int m, nnz, *csrRowPtr = NULL, *csrColInd = NULL;
	if (strstr(argv[2], ".mtx"))
		mtx2csr(argv[2], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[2], ".graph"))
		graph2csr(argv[2], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[2], ".gr"))
		gr2csr(argv[2], m, nnz, csrRowPtr, csrColInd);
	else
		{ printf("Unrecognizable input file format\n"); exit(0); }
	int *coloring, correct;
	coloring = (int *)calloc(m, sizeof(int));
	correct = 1;
	color(m, nnz, csrRowPtr, csrColInd, coloring);
	verify(m, nnz, csrRowPtr, csrColInd, coloring, &correct);
	if (correct)
		printf("correct\n");
	else
		printf("incorrect\n");
	write_solution("coloring.txt", m, coloring);
	return 0;
}
