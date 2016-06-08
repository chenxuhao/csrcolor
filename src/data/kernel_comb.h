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
#define	MAXCOLOR 128
typedef cub::BlockScan<int, BLKSIZE> BlockScan;

__global__ void initialize(int *coloring, int m) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < m) {
		coloring[id] = MAXCOLOR;
	}   
}

__device__ __forceinline__ void assignColor(unsigned *forbiddenColors, int *coloring, int vertex) {
	int vertex_color;
	for (vertex_color = 0; vertex_color < MAXCOLOR/32; vertex_color++) {
		int pos = __ffs(forbiddenColors[vertex_color]);
		if(pos) {
			coloring[vertex] = vertex_color * 32 + pos - 1;
			break;
		}
	}
	assert(vertex_color < MAXCOLOR);
}

__global__ void firstFit(int m, const int* __restrict__ csrRowPtr, const int* __restrict__ csrColInd, Worklist2 inwl, int *coloring) {
//__global__ void firstFit(int m, int* csrRowPtr, int* csrColInd, Worklist2 inwl, int *coloring) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned forbiddenColors[MAXCOLOR/32+1];
	int id = tid;
	//int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	//for (int id = tid; total_inputs > 0; id += blockDim.x * gridDim.x, total_inputs--) {
		int vertex;
		if (inwl.pop_id(id, vertex)) {
			//int row_begin = csrRowPtr[vertex];
			//int row_end = csrRowPtr[vertex + 1]; 
			int row_begin = __ldg(csrRowPtr + vertex);
			int row_end= __ldg(csrRowPtr + vertex + 1);
			for (int j = 0; j < MAXCOLOR/32; j++)
				forbiddenColors[j] = 0xffffffff;
			for (int offset = row_begin; offset < row_end; offset ++) {
				//int neighbor = csrColInd[offset];
				int neighbor = __ldg(csrColInd + offset);
				int color = coloring[neighbor];
				forbiddenColors[color / 32] &= ~(1 << (color % 32));
			}
			assignColor(forbiddenColors, coloring, vertex);
		}
	//}
}

__device__ __forceinline__ void conflictDetect1(int src, int dst, int *coloring, bool &is_conflict) {
    int color_s = coloring[src];
	int color_d = coloring[dst];
	if (color_s == color_d && src < dst) {
		is_conflict = 1;
		coloring[src] = MAXCOLOR;
	}
}

__device__ __forceinline__ bool conflictDetect2(int src, int dst, int *coloring, int *degree, bool &is_conflict) {
	if (coloring[src] == coloring[dst]) {
		bool is_victim;
		if (degree[src] == degree[dst])
			is_victim = (src < dst) ? true : false;
		else is_victim = (degree[src] < degree[dst]) ? true : false;
		if (is_victim) {
			is_conflict = 1;
			coloring[src] = MAXCOLOR;
		}
	}
}

__global__ void conflictResolve(int m, const int* __restrict__ csrRowPtr, const int* __restrict__ csrColInd, Worklist2 inwl, Worklist2 outwl, int * degree, int *coloring) {
//__global__ void conflictResolve(int m, int* csrRowPtr, int* csrColInd, Worklist2 inwl, Worklist2 outwl, int * degree, int *coloring) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	bool is_conflict = 0;
	int vertex;
	if (inwl.pop_id(id, vertex)) {
		//int row_begin = csrRowPtr[vertex];
		//int row_end = csrRowPtr[vertex + 1]; 
		int row_begin = __ldg(csrRowPtr + vertex);
		int row_end= __ldg(csrRowPtr + vertex + 1);
		for (int offset = row_begin; offset < row_end; offset ++) {
			//int neighbor = csrColInd[offset];
			int neighbor = __ldg(csrColInd + offset);
			//conflictDetect1(vertex, neighbor, coloring, is_conflict);
			conflictDetect2(vertex, neighbor, coloring, degree, is_conflict);
			if(is_conflict) break;
		}
	}   
	outwl.push_1item<BlockScan>((int)is_conflict, vertex, BLKSIZE);
}

void color(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int num_SMs) {
	double starttime, endtime;
	double runtime[ITERATIONS];
	int colors[ITERATIONS];
	int iterations[ITERATIONS];
	int *d_csrRowPtr, *d_csrColInd, *d_coloring, *d_degree;
	printf("Graph coloring data-driven Combination A version\n");
	int *degree = (int *)malloc(m * sizeof(int));
	for(int i = 0; i < m; i ++) {
		degree[i] = csrRowPtr[i + 1] - csrRowPtr[i];
	}
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrRowPtr, (m + 1) * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrColInd, nnz * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_degree, m * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_coloring, m * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrRowPtr, csrRowPtr, (m + 1) * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrColInd, csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaDeviceSynchronize());
	int device = 0;
	int deviceCount = 0;
	CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	int nSM = deviceProp.multiProcessorCount;
	//int nSM = num_SMs;
	fprintf(stdout, "Found %d devices, using device %d (%s), compute capability %d.%d, cores %d*%d.\n",
		deviceCount, device, deviceProp.name, deviceProp.major, deviceProp.minor, nSM, ConvertSMVer2Cores(deviceProp.major, deviceProp.minor));
	
	const size_t max_blocks = maximum_residency(firstFit, BLKSIZE, 0); 
	printf("max_blocks=%d\n", max_blocks);

	for (int i = 0; i < ITERATIONS; i++) {
		Worklist2 inwl(m), outwl(m);
		Worklist2 *inwlptr = &inwl, *outwlptr = &outwl;
		CUDA_SAFE_CALL(cudaMemcpy(inwl.dindex, &m, sizeof(int), cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(d_degree, degree, m * sizeof(int), cudaMemcpyHostToDevice));
		initialize <<<((m - 1) / BLKSIZE + 1), BLKSIZE>>> (d_coloring, m);
		CUDA_SAFE_CALL(cudaDeviceSynchronize());
		iterations[i] = 0;

		starttime = rtclock();
		int nitems = m;
		thrust::sequence(thrust::device, inwl.dwl, inwl.dwl + m);
		while (nitems > 0) {
			iterations[i] ++;
			//printf("in_nitems[%d]=%d\n", iteration, nitems);
			int nblocks = (nitems - 1) / BLKSIZE + 1;
			int nblocks_1 = nSM * max_blocks;
			if(nblocks < nblocks_1) nblocks_1 = nblocks;
			firstFit<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, *inwlptr, d_coloring);
			conflictResolve<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, *inwlptr, *outwlptr, d_degree, d_coloring);
			nitems = outwlptr->nitems();
			Worklist2 * tmp = inwlptr;
			inwlptr = outwlptr;
			outwlptr = tmp;
			outwlptr->reset();
		}
		CUDA_SAFE_CALL(cudaDeviceSynchronize());
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
