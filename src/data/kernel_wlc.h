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

__global__ void initialize(int *coloring, int m) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < m) {
		coloring[id] = MAXCOLOR;
	}   
}

__global__ void FirstFit(int m, int *csrRowPtr, int *csrColInd, Worklist2 inwl, int *coloring) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	bool forbiddenColors[MAXCOLOR+1];
	int vertex;
	if (inwl.pop_id(id, vertex)) {
		for (int j = 0; j < MAXCOLOR; j++)
			forbiddenColors[j] = false;
		int row_begin = csrRowPtr[vertex];
		int row_end = csrRowPtr[vertex + 1];
		for (int offset = row_begin; offset < row_end; offset ++) {
			int neighbor = csrColInd[offset];
			int color = coloring[neighbor];
			forbiddenColors[color] = true;
		}
		int vertex_color;
		for (vertex_color = 0; vertex_color < MAXCOLOR; vertex_color ++) {
			if (!forbiddenColors[vertex_color]) {
				coloring[vertex] = vertex_color;
				break;
			}
		}   
		assert(vertex_color < MAXCOLOR);
	}
}

__device__ __forceinline__ bool conflictDetect(int src, int dst, int *coloring, bool &is_conflict) {
	//int color_s = coloring[src];
	int color_s = cub::ThreadLoad<cub::LOAD_CG>(coloring + src);
	//int color_d = coloring[dst];
	int color_d = cub::ThreadLoad<cub::LOAD_CG>(coloring + dst);
	if (color_s == color_d && src < dst) {
		is_conflict = 1;
		//coloring[src] = MAXCOLOR;
		cub::ThreadStore<cub::STORE_CG>(coloring + src, MAXCOLOR);
		return true;
	}
	return is_conflict;
}

__device__ __forceinline__ unsigned LaneId() {
	unsigned ret;
	asm("mov.u32 %0, %laneid;" : "=r"(ret));
	return ret;
}

__device__ __forceinline__ void resolveByCta(int tid, int m, int vertex, int size, int *row_offsets, int *column_indices, int *degree, Worklist2 &inwl, int *coloring, bool *is_conflict) {
	__shared__ int owner;
	__shared__ int sh_vertex;
	owner = -1;
	int id = tid;
	int owner_tid;
	int my_size = size;
	//int vertex;
	//if(inwl.pop_id(id, vertex)) {
		//size = row_end - row_begin;
		//size = degree[vertex];
	//}
	while(true) {
		if(my_size >= BLKSIZE)
			owner = threadIdx.x;
		__syncthreads();
		owner_tid = owner;
		if(owner_tid == -1)
			break;
		__syncthreads();
		if(owner == threadIdx.x) {
			sh_vertex = vertex;
			// mark this vertex as processed already
			//inwl.dwl[id] = -1;
			cub::ThreadStore<cub::STORE_CG>(inwl.dwl + id, -1);
			owner = -1;
			my_size = 0;
		}
		__syncthreads();
		int row_begin = row_offsets[sh_vertex];
		//row_begin = __ldg(row_offsets + sh_vertex);
		//row_end = row_offsets[sh_vertex + 1];
		//row_end = __ldg(row_offsets + sh_vertex + 1);
		//int neighbor_size = row_end - row_begin;
		int neighbor_size = degree[sh_vertex];
		int num = ((neighbor_size - 1) / BLKSIZE + 1) * BLKSIZE;
		for(int i = threadIdx.x; i < num; i += BLKSIZE) {
			if(i < neighbor_size) {
				int edge = row_begin + i;
				int dst = column_indices[edge];
				//int dst = __ldg(column_indices + edge);
				if(conflictDetect(sh_vertex, dst, coloring, is_conflict[owner_tid]))
					break;
			}
		}
	}
}

#define WARP_SIZE 32
#define LOG_WARP_SIZE 5
#define NUM_WARPS (BLKSIZE / WARP_SIZE)
__device__ __forceinline__ void resolveByWarp(int id, int m, int *row_offsets, int *column_indices, int *degree, Worklist2 &inwl, int *coloring, bool *is_conflict) {
	unsigned warp_id = threadIdx.x >> LOG_WARP_SIZE;
	unsigned lane_id = LaneId();
	__shared__ int owner[NUM_WARPS];
	__shared__ int sh_vertex[NUM_WARPS];
	owner[warp_id] = -1;
	int size = 0;
	int vertex;
	if(inwl.pop_id(id, vertex)) {
		if (vertex != -1) {
			//size = row_end - row_begin;
			size = degree[vertex];
		}
	}
	while(__any(size) >= WARP_SIZE) {
		if(size >= WARP_SIZE)
			owner[warp_id] = lane_id;
		if(owner[warp_id] == lane_id) {
			sh_vertex[warp_id] = vertex;
			// mark this vertex as processed already
			//inwl.dwl[id] = -1;
			cub::ThreadStore<cub::STORE_CG>(inwl.dwl + id, -1);
			owner[warp_id] = -1;
			size = 0;
		}
		int winner = sh_vertex[warp_id];
		int winner_tid = warp_id * WARP_SIZE + owner[warp_id];
		int row_begin = row_offsets[winner];
		//row_begin = __ldg(row_offsets + sh_vertex[warp_id]);
		//row_end = row_offsets[sh_vertex[warp_id] + 1];
		//row_end = __ldg(row_offsets + sh_vertex[warp_id] + 1);
		//int neighbor_size = row_end - row_begin;
		int neighbor_size = degree[winner];
		int num = ((neighbor_size + WARP_SIZE - 1) / WARP_SIZE) * WARP_SIZE;
		for(int i = lane_id; i < num; i+= WARP_SIZE) {
			if(i < neighbor_size) {
				int edge = row_begin + i;
				int dst = column_indices[edge];
				//int dst = __ldg(column_indices + edge);
				if(conflictDetect(winner, dst, coloring, is_conflict[winner_tid])) break;
			}
		}
	}
}

__global__ void conflictResolve(int m, int *row_offsets, int *column_indices, int *degree, Worklist2 inwl, Worklist2 outwl, int *coloring) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tx = threadIdx.x;
	__shared__ bool is_conflict[BLKSIZE];
	is_conflict[tx] = 0;
	__syncthreads();
	int id = tid;
	int old_vertex;
	int my_degree = 0;
	if(inwl.pop_id(id, old_vertex)) {
		my_degree = degree[old_vertex];
	}

	resolveByCta(id, m, old_vertex, my_degree, row_offsets, column_indices, degree, inwl, coloring, is_conflict);
	//resolveByWarp(id, m, row_offsets, column_indices, degree, inwl, coloring, is_conflict);

	typedef cub::BlockScan<int, BLKSIZE> BlockScan;
	__shared__ BlockScan::TempStorage temp_storage;
	__shared__ int gather_offsets[SCRATCHSIZE];
	__shared__ int src[BLKSIZE];
	__shared__ short srcIndex[BLKSIZE];
	gather_offsets[tx] = 0;
	src[tx] = 0;
	srcIndex[tx] = 0;

	//int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	//for (int id = tid; total_inputs > 0; id += blockDim.x * gridDim.x, total_inputs--) {
	int vertex;
	int neighbor_size = 0;
	int row_begin = 0;
	//int row_end = 0;
	int scratch_offset = 0;
	int total_edges = 0;
	if (inwl.pop_id(id, vertex)) {
		if (vertex != -1) {
			row_begin = row_offsets[vertex];
			//row_begin = __ldg(row_offsets + vertex);
			//row_end = row_offsets[vertex + 1];
			//row_end = __ldg(row_offsets + vertex + 1);
			//neighbor_size = row_end - row_begin;
			neighbor_size = degree[vertex];
		}
	}
	BlockScan(temp_storage).ExclusiveSum(neighbor_size, scratch_offset, total_edges);
	int done = 0;
	int neighborsdone = 0;

	while (total_edges > 0) {
		__syncthreads();
		int i;
		for (i = 0; (neighborsdone + i) < neighbor_size && (scratch_offset + i - done) < SCRATCHSIZE; i++) {
			int ii = scratch_offset + i - done;
			gather_offsets[ii] = row_begin + neighborsdone + i;
			src[ii] = vertex;
			srcIndex[ii] = tx;
		}
		neighborsdone += i;
		scratch_offset += i;
		__syncthreads();

		if (tx < total_edges) {
			int edge = gather_offsets[tx];
			int dst = column_indices[edge];
			//int dst = __ldg(column_indices + edge);
			int index = srcIndex[tx];
			int srcsrc = src[tx];
			conflictDetect(srcsrc, dst, coloring, is_conflict[index]);
		}
		total_edges -= BLKSIZE;
		done += BLKSIZE;
	}
	__syncthreads();
	outwl.push_1item<BlockScan>((int)is_conflict[tx], old_vertex, BLKSIZE);
	//}
}

void color(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int num_SMs) {
	double starttime, endtime;
	double runtime[ITERATIONS];
	int colors[ITERATIONS];
	int iterations[ITERATIONS];
	int *d_csrRowPtr, *d_csrColInd, *d_coloring, *d_degree;
	printf("Graph coloring data-driven Load Balancing version\n");

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
	CUDA_SAFE_CALL(cudaMemcpy(d_degree, degree, m * sizeof(int), cudaMemcpyHostToDevice));
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
		while (nitems > 0) {
			iterations[i] ++;
			//printf("in_nitems[%d]=%d\n", iteration, nitems);
			int nblocks = (nitems - 1) / BLKSIZE + 1;
			//printf("Iteration %d (nitems=%d, nblocks=%d, blocksize=%d).\n", iteration, nitems, nblocks, BLKSIZE);
			FirstFit<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, *inwlptr, d_coloring);
			conflictResolve<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, d_degree, *inwlptr, *outwlptr, d_coloring);
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
