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
#include "worklistc.h"
#define	SCRATCHSIZE BLKSIZE
#define	MAXCOLOR 128 // assume graph can be colored with less than 128 colors
//#define TEXTURE

typedef cub::BlockScan<int, BLKSIZE> BlockScan;
#ifdef TEXTURE
texture <int, 1, cudaReadModeElementType> rowPtr;
texture <int, 1, cudaReadModeElementType> colInd;
#endif

__global__ void initialize(int *coloring, int m) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < m) {
		coloring[id] = MAXCOLOR;
	}   
}

#ifdef TEXTURE
__global__ void FirstFit(int m, Worklist2 inwl, int *coloring) {
#else
__global__ void FirstFit(int m, const int * __restrict__  csrRowPtr, const int * __restrict__ csrColInd, Worklist2 inwl, int *coloring) {
#endif
	int id = blockIdx.x * blockDim.x + threadIdx.x;	
	bool forbiddenColors[MAXCOLOR+1];
	int vertex;
	if (inwl.pop_id(id, vertex)) {
#ifdef TEXTURE
		int row_begin = tex1Dfetch(rowPtr, vertex);
		int row_end = tex1Dfetch(rowPtr, vertex + 1);
#else
		int row_begin = __ldg(csrRowPtr + vertex);
		int row_end = __ldg(csrRowPtr + vertex + 1);
#endif
		for (int i = 0; i < MAXCOLOR; i ++)
			forbiddenColors[i] = false;
		for (int offset = row_begin; offset < row_end; offset ++) {
#ifdef TEXTURE
			int neighbor = tex1Dfetch(colInd, offset);
#else
			int neighbor = __ldg(csrColInd + offset);
#endif
			int color = coloring[neighbor];
			forbiddenColors[color] = true;
		}
		int vertex_color;
		for (vertex_color = 0; vertex_color < MAXCOLOR; vertex_color++) {
			if (!forbiddenColors[vertex_color]) {
				coloring[vertex] = vertex_color;
				break;
			}
		}
		assert(vertex_color < MAXCOLOR);
	}
}

#ifdef TEXTURE
__global__ void conflictResolve(int m, Worklist2 inwl, Worklist2 outwl, int *coloring) {
#else
__global__ void conflictResolve(int m, const int * __restrict__  csrRowPtr, const int * __restrict__  csrColInd, Worklist2 inwl, Worklist2 outwl, int *coloring) {
#endif
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int conflicted = 0;
	int vertex;
	if (inwl.pop_id(id, vertex)) {
#ifdef TEXTURE
		int row_begin = tex1Dfetch(rowPtr, vertex);
		int row_end= tex1Dfetch(rowPtr, vertex + 1);
#else
		int row_begin = __ldg(csrRowPtr + vertex);
		int row_end= __ldg(csrRowPtr + vertex + 1);
#endif
		for (int offset = row_begin; offset < row_end; offset ++) {
#ifdef TEXTURE
			int neighbor = tex1Dfetch(colInd, offset);
#else
			int neighbor = __ldg(csrColInd + offset);
#endif
			if (coloring[vertex] == coloring[neighbor] && vertex < neighbor) {
				conflicted = 1;
				coloring[vertex] = MAXCOLOR;
				break;
			}
		}
	}
	outwl.push_1item<BlockScan>(conflicted, vertex, BLKSIZE);
}

void color(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int num_SMs) {
	double starttime, endtime;
	double runtime[ITERATIONS];
	int colors[ITERATIONS];
	int iterations[ITERATIONS];
	int *d_csrRowPtr, *d_csrColInd, *d_coloring;
	printf("Graph coloring data-driven LDG version\n");
	
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrRowPtr, (m + 1) * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrColInd, nnz * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_coloring, m * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrRowPtr, csrRowPtr, (m + 1) * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrColInd, csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice));
#ifdef TEXTURE
	CUDA_SAFE_CALL(cudaBindTexture(0, rowPtr, csrRowPtr, (m + 1) * sizeof(int)));
	CUDA_SAFE_CALL(cudaBindTexture(0, colInd, csrColInd, (nnz + 1) * sizeof(int)));
#endif
	CUDA_SAFE_CALL(cudaDeviceSynchronize());
	int device = 0;
	int deviceCount = 0;
	CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, device);
	int nSM = deviceProp.multiProcessorCount;
	fprintf(stdout, "Found %d devices, using device %d (%s), compute capability %d.%d, cores %d*%d.\n", 
			deviceCount, device, deviceProp.name, deviceProp.major, deviceProp.minor, nSM, ConvertSMVer2Cores(deviceProp.major, deviceProp.minor));

	for (int i = 0; i < ITERATIONS; i++) {
		Worklist2 inwl(m), outwl(m);
		Worklist2 *inwlptr = &inwl, *outwlptr = &outwl;
		CUDA_SAFE_CALL(cudaMemcpy(inwl.dindex, &m, sizeof(int), cudaMemcpyHostToDevice));
		initialize <<<((m - 1) / BLKSIZE + 1), BLKSIZE>>> (d_coloring, m);
		CUDA_SAFE_CALL(cudaDeviceSynchronize());
		iterations[i] = 0;

		starttime = rtclock();
		int nitems = m;
		thrust::sequence(thrust::device, inwl.dwl, inwl.dwl + m);
		int iteration = 0;
		while (nitems > 0) {
			iterations[i] ++;
			int nblocks = (nitems - 1) / BLKSIZE + 1;
#ifdef TEXTURE
			FirstFit<<<nblocks, BLKSIZE>>>(m, *inwlptr, d_coloring);
			conflictResolve<<<nblocks, BLKSIZE>>>(m, *inwlptr, *outwlptr, d_coloring);
#else
			FirstFit<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, *inwlptr, d_coloring);
			conflictResolve<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, *inwlptr, *outwlptr, d_coloring);
#endif
			nitems = outwlptr->nitems();
			Worklist2 * tmp = inwlptr;
			inwlptr = outwlptr;
			outwlptr = tmp;
			outwlptr->reset();
		}
		cudaDeviceSynchronize();
		endtime = rtclock();
		runtime[i] = 1000.0f * (endtime - starttime);
		colors[i] = thrust::reduce(thrust::device, d_coloring, d_coloring + m, 0, thrust::maximum<int>()) + 1;
	}
	CUDA_SAFE_CALL(cudaMemcpy(coloring, d_coloring, m * sizeof(int), cudaMemcpyDeviceToHost));
	double total_time = 0.0;
	int total_colors = 0;
	int total_iterations = 0;
	for (int i = 0; i < ITERATIONS; i++) {
		total_time += runtime[i];
		total_colors += colors[i];
		total_iterations += iterations[i];
		printf("[%d %.2f %d] ", colors[i], runtime[i], iterations[i]);
	}   
	double avg_time = (double)total_time / ITERATIONS;
	double avg_colors = (double)total_colors / ITERATIONS;
	double avg_iterations = (double)total_iterations / ITERATIONS;
	printf("\navg_time %f ms, avg_colors %.2f avg_iterations %.2f\n", avg_time, avg_colors, avg_iterations);
	CUDA_SAFE_CALL(cudaFree(d_csrRowPtr));
	CUDA_SAFE_CALL(cudaFree(d_csrColInd));
	CUDA_SAFE_CALL(cudaFree(d_coloring));
}
