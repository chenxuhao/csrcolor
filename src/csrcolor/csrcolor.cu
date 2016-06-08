// Copyright 2016, National University of Defense Technology
// Authors: Xuhao Chen <cxh@illinois.edu> and Pingfan Li <lipingfan@163.com>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "cusparse.h"
#include "cuda.h"
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include "cutil_subset.h"
#include "common.h"
#include <vector>
#include <set>
using namespace std;

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
	//printf("%c %s %d %d\n", c, sp, m, nnz);
	printf("num_vertices %d  num_edges %d\n", m, nnz);
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

void write_solution(char *fname, int m, int *coloring) {
	FILE *fp = fopen(fname, "w");
	int i;
	for (i = 0; i < m; i++) {
		fprintf(fp, "%d\n", coloring[i]);
	}
	fclose(fp);
}

void verify(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int *correct) {
	int i, j, neighbors, start, neighbor_j;
	for (i = 0; i < m; i++) {
		start = csrRowPtr[i];
		neighbors = csrRowPtr[i + 1] - start;
		for (j = 0; j < neighbors; j++) {
			neighbor_j = csrColInd[start + j];
			if (coloring[i] == coloring[neighbor_j] && i != neighbor_j) {
				*correct = 0;
				printf("coloring[%d] = coloring[%d] = %d\n", i, neighbor_j, coloring[i]);
			}
			break;
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
	printf("num_vertices %d, num_edges %d\n", m, nnz);
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


int main(int argc, char *argv[]) {
	if (argc != 2) {
		printf("Usage: %s <graph>\n", argv[0]);
		exit(1);
	}
	int m, nnz, *csrRowPtr = NULL, *csrColInd = NULL;
	if (strstr(argv[1], ".mtx"))
		mtx2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[1], ".graph"))
		graph2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else if (strstr(argv[1], ".gr"))
		gr2csr(argv[1], m, nnz, csrRowPtr, csrColInd);
	else
		{ printf("Unrecognizable input file format\n"); exit(0); }
	if (csrRowPtr == NULL)
		printf("csrRowPtr is NULL\n");
	if (csrColInd == NULL)
		printf("csrColInd is NULL\n");
	double t1, t2, t3, t4, t5, t6;
	int *d_csrRowPtr, *d_csrColInd;
	float *d_csrVal;
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrRowPtr, (m + 1) * sizeof(int))); 
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrColInd, nnz * sizeof(int))); 
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrVal, nnz * sizeof(float))); 
	int ncolors = 0, *coloring;
	int *d_coloring, *d_reordering;
	float fraction = 1.0;
	coloring = (int *)calloc(m, sizeof(int));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_coloring, m * sizeof(int)));	
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_reordering, m * sizeof(int))); 
	CUDA_SAFE_CALL(cudaMemset(d_reordering, 0, m * sizeof(int))); 
	CUDA_SAFE_CALL(cudaDeviceSynchronize());
	t1 = rtclock();
	CUDA_SAFE_CALL(cudaMemcpy(d_csrRowPtr, csrRowPtr, (m + 1) * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrColInd, csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	t2 = rtclock();
	//printf("time of init:%f ms\n", 1000.0f * (t2 - t1));	

	int device = 0;
	int deviceCount = 0;
	CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
	cudaDeviceProp deviceProp;
	CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, device));
	int nSM = deviceProp.multiProcessorCount;
	fprintf(stdout, "Found %d devices, using device %d (%s), compute capability %d.%d, cores %d*%d.\n", 
			deviceCount, device, deviceProp.name, deviceProp.major, deviceProp.minor, nSM, ConvertSMVer2Cores(deviceProp.major, deviceProp.minor));

	cusparseStatus_t status;
	cusparseHandle_t handle;
	status = cusparseCreate(&handle);
	if (status != CUSPARSE_STATUS_SUCCESS) {
		printf("error!");
		exit(1);
	}
	cusparseMatDescr_t descr;
	status = cusparseCreateMatDescr(&descr);
	if (status != CUSPARSE_STATUS_SUCCESS) {
		printf("error!");
		exit(1);
	}
	cusparseColorInfo_t info;
	status = cusparseCreateColorInfo(&info);
	if (status != CUSPARSE_STATUS_SUCCESS) {
		printf("error!");
		exit(1);
	}	
	double runtime[10];
	int colors[10];
	for (int i = 0; i < 10; i++) {
		t5 = rtclock();
		status = cusparseScsrcolor(handle, m, nnz, descr, d_csrVal, d_csrRowPtr, d_csrColInd, &fraction, &ncolors, d_coloring, d_reordering, info);
		t6 = rtclock();
		runtime[i] = 1000.0f * (t6 - t5);
		colors[i] = 1 + thrust::reduce(thrust::device, d_coloring, d_coloring + m, 0, thrust::maximum<int>());
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
	switch (status) {
		case CUSPARSE_STATUS_SUCCESS:
			//printf("success\n");
			break;
		case CUSPARSE_STATUS_NOT_INITIALIZED:
			printf("not initialed\n");
		case CUSPARSE_STATUS_ALLOC_FAILED:
			printf("alloc failed\n");
			break;
		case CUSPARSE_STATUS_INVALID_VALUE:
			printf("invalid value\n");
			break;
		case CUSPARSE_STATUS_ARCH_MISMATCH:
			printf("mismatch\n");
			break;
		case CUSPARSE_STATUS_INTERNAL_ERROR:
			printf("internal error\n");
			break;
		case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
			printf("not supported\n");
			break;
		default:
			printf("unknown error\n");
			break;
	};
	t3 = rtclock();
	CUDA_SAFE_CALL(cudaMemcpy(coloring, d_coloring, m * sizeof(int), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaDeviceSynchronize());
	t4 = rtclock();
	//printf("time of copy back:%f ms\n", 1000.0f * (t4 - t3));	
	write_solution("color.txt", m, coloring);
	int correct = 1;
	verify(m, nnz, csrRowPtr, csrColInd, coloring, &correct);
	if (correct)
		printf("correct.\n");
	else
		printf("incorrect.\n");
	return 0;
}
