// Copyright 2016, National University of Defense Technology
// Authors: Xuhao Chen <cxh@illinois.edu> and Pingfan Li <lipingfan@163.com>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
using namespace std;
#include "graph_io.h"
/*
// transfer *.graph file to CSR format
void graph2csr(char *graph, int &m, int &nnz, int *&csrRowPtr, int *&csrColInd) {
	printf("Reading .graph input file %s\n", graph);
	std::ifstream cfile;
	cfile.open(graph);
	std::string str;
	getline(cfile, str);
	sscanf(str.c_str(), "%d %d", &m, &nnz);
	printf("num_vertices %d  num_edges %d\n", m, nnz);
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
	printf("num_vertices %d  num_edges %d\n", m, nnz);
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
	printf("Reading .mtx input file %s\n", mtx);
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
	printf("num_vertices %d  num_edges %d\n", m, nnz);
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

void write_solution(char *fname, int *coloring, int n) {
	int i;
	FILE *fp;
	fp = fopen(fname, "w");
	for (i = 0; i < n; i++) {
		fprintf(fp, "%d\n", coloring[i]);
	}
	fclose(fp);
}

void verify(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int *correct) {
	int i, offset, neighbor_j;
	for (i = 0; i < m; i++) {
		for (offset = csrRowPtr[i]; offset < csrRowPtr[i + 1]; offset++) {
			neighbor_j = csrColInd[offset];
			if (coloring[i] == coloring[neighbor_j] && neighbor_j != i) {
				*correct = 0;
				printf("coloring[%d] = coloring[%d] = %d\n", i, neighbor_j, coloring[i]);
				return;
			}
		}	
	}
}
//*/

double rtclock() {
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

#define MAXCOLOR 128
void FirstFit(int m, int nnz, int *csrRowPtr, int *csrColInd, int *ncolors, int *coloring) {
	int max_color = 1;
	int vertex;
	int forbiddenColors[MAXCOLOR+1];
	for (int i = 0; i < MAXCOLOR; i ++)
		forbiddenColors[i] = -1;
	for (vertex = 0; vertex < m; vertex++) {
		int row_begin = csrRowPtr[vertex];
		int row_end = csrRowPtr[vertex + 1];
		for (int offset = row_begin; offset < row_end; offset++) {
			int neighbor = csrColInd[offset];
			forbiddenColors[coloring[neighbor]] = vertex;
		}
		int vertex_color = 1;
		while (vertex_color < max_color && forbiddenColors[vertex_color] == vertex)
			vertex_color++;
		if (vertex_color == max_color)
			max_color++;
		assert(vertex_color < MAXCOLOR);
		coloring[vertex] = vertex_color;
	}
	*ncolors = max_color - 1;
}

int main(int argc, char *argv[]) {
	int iteration = 0;
	int m, nnz, *csrRowPtr = NULL, *csrColInd = NULL;
	if (strstr(argv[1], ".mtx"))
		mtx2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[1], ".graph"))
		graph2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[1], ".gr"))
		gr2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else
		{ printf("Unrecognizable input file format\n"); exit(0); }
	int ncolors, *coloring, correct;
	ncolors = 0;
	coloring = (int *)calloc(m, sizeof(int));
	correct = 1;
	double starttime, endtime;
	double runtime[10];
	int colors[10];
	for (int i = 0; i < 10; i++) {
		memset(coloring, 0, m * sizeof(int));	
		starttime = rtclock();
		FirstFit(m, nnz, csrRowPtr, csrColInd, &ncolors, coloring);
		endtime = rtclock();
		runtime[i] = (1000.0f) * (endtime - starttime);
		colors[i] = ncolors;
	}
	double total_time = 0;
	int total_colors = 0;
	double avg_time;
	double avg_colors;
	for (int i = 0; i < 10; i++) {
		printf("[%.2f %d] ", runtime[i], colors[i]);
		total_time += runtime[i];
		total_colors += colors[i];
	}
	printf("\navg_time %f ms, avg_colors %.2f\n", total_time / 10, (double)total_colors / 10);
	write_solution("color.txt", coloring, m);
	verify(m, nnz, csrRowPtr, csrColInd, coloring, &correct);
	if (correct)
		printf("correct.\n");
	else
		printf("incorrect.\n");
	return 0;
}
