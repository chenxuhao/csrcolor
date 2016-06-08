// Copyright 2016, National University of Defense Technology
// Authors: Xuhao Chen <cxh@illinois.edu> and Pingfan Li <lipingfan@163.com>
#include <cub/cub.cuh>
#include "gbar.cuh"
#include "cuda_launch_config.hpp"
#include "cutil_subset.h"
#include "common.h"
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/count.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <b40c_test_util.h>
#include <b40c/graph/builder/dimacs.cuh>
#include <b40c/graph/color/csr_problem.cuh>
#include <b40c/graph/csr_graph.cuh>
//#include <b40c/graph/color/enactor_hybrid.cuh>
#include <b40c/graph/color/enactor_two_phase.cuh>
#define INIT_VAL -1
using namespace b40c;
using namespace graph;

void verify_color(unsigned *dist, int m, int *csrRowPtr, int *csrColInd, unsigned *nerr) {
	for (int nn = 0; nn < m; nn ++) {
		int neighbor_offset = csrRowPtr[nn];
		int neighbor_size = csrRowPtr[nn + 1] - neighbor_offset;
		for (unsigned ii = 0; ii < neighbor_size; ++ii) {
			int v = csrColInd[neighbor_offset + ii];
			unsigned wt = 1;
			if (wt > 0 && dist[nn] + wt < dist[v]) {
				//printf("%d %d %d %d\n", nn, v, dist[nn], dist[v]);
				++*nerr;
			}
		}
	}	
}

void write_solution(const char *fname, int m, unsigned *h_dist) {
	//unsigned *h_dist;
	//h_dist = (unsigned *) malloc(m * sizeof(unsigned));
	assert(h_dist != NULL);
	//CUDA_SAFE_CALL(cudaMemcpy(h_dist, dist, m * sizeof(foru), cudaMemcpyDeviceToHost));
	printf("Writing solution to %s\n", fname);
	FILE *f = fopen(fname, "w");
	fprintf(f, "Computed solution (source dist): [");
	for(int node = 0; node < m; node++) {
		fprintf(f, "%d:%d\n ", node, h_dist[node]);
	}
	fprintf(f, "]");
	free(h_dist);
}

void color_ldb(int m, int nnz, int *csrRowPtr, int *csrColInd, int *ncolors, int *coloring, int num_SMs) {
	printf("Graph coloring data-driven load-balance version\n");
	typedef int VertexId;
	typedef unsigned Value;
	typedef int SizeT;
	int *d_csrRowPtr, *d_csrColInd;
	int *d_coloring;
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrRowPtr, (m + 1) * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrColInd, nnz * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_coloring, m * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrRowPtr, csrRowPtr, (m + 1) * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrColInd, csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemset(d_coloring, INIT_VAL, m  * sizeof(int)));

	graph::CsrGraph<VertexId, Value, SizeT> csr_graph;
	csr_graph.FromScratch<true>(m, nnz);
	CUDA_SAFE_CALL(cudaMemcpy(csr_graph.row_offsets, d_csrRowPtr, sizeof(SizeT) * (m + 1), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(csr_graph.column_indices, d_csrColInd, sizeof(VertexId) * nnz, cudaMemcpyDeviceToHost));

	typedef color::CsrProblem<VertexId, SizeT, false> CsrProblem;
	color::EnactorTwoPhase<false> two_phase(false);
	//color::EnactorHybrid<false> hybrid(false);
	CsrProblem csr_problem;
	if (csr_problem.FromHostProblem(false, csr_graph.nodes, csr_graph.edges, csr_graph.column_indices, csr_graph.row_offsets, 1)) exit(1);
	cudaError_t	retval = cudaSuccess;
	double runtime[ITERATIONS];
	int colors[ITERATIONS];
	double starttime, endtime;
	for (int i = 0; i < ITERATIONS; i++) {
		starttime = rtclock();
		//if (retval = csr_problem.Reset(hybrid.GetFrontierType(), 1.3))
		if (retval = csr_problem.Reset(two_phase.GetFrontierType(), 1.3))
			return;
		//if (retval = hybrid.EnactSearch(csr_problem, 0)) {
		if (retval = two_phase.EnactIterativeSearch(csr_problem, 0)) {
			if (retval && (retval != cudaErrorInvalidDeviceFunction)) {
				exit(1);
			}
		}
		CUDA_SAFE_CALL(cudaDeviceSynchronize());
		endtime = rtclock();
		runtime[i] = 1000.0f * (endtime - starttime);
		//colors[i] = thrust::reduce(thrust::device, d_coloring, d_coloring + m, 0, thrust::maximum<int>());
	}
	unsigned *h_dist;
	h_dist = (unsigned *) malloc(m * sizeof(unsigned));
	assert(h_dist != NULL);
	if (csr_problem.ExtractResults((int *) h_dist)) exit(1);
	for(int i = 0; i < m; i++)
		if((signed) h_dist[i] == -1)
			h_dist[i] = 1000000000;
	printf("Done!\n");
	unsigned nerr = 0;
	printf("verifying.\n");
	verify_color(h_dist, m, csrRowPtr, csrColInd, &nerr);
	printf("\tno of errors = %d.\n", nerr);
	write_solution("color-output.txt", m, h_dist);
	exit(0);
	/*
	cudaMemcpy(coloring, d_coloring, m * sizeof(int), cudaMemcpyDeviceToHost);
	*ncolors = colors[ITERATIONS - 1];
	double totaltime = 0.0;
	int totalcolors = 0;
	for (int i = 0; i < ITERATIONS; i++) {
		totaltime += runtime[i];
		totalcolors += colors[i];
		printf("[%d %f] ", colors[i], runtime[i]);
	}
	double avgtime = (double)totaltime / ITERATIONS;
	double avgcolors = (double)totalcolors / ITERATIONS;
	printf("\navgtime=%f ms, avgcolors = %f\n", avgtime, avgcolors);
	*/
}
