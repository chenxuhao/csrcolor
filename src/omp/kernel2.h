// Copyright 2016, National University of Defense Technology
// Authors: Xuhao Chen <cxh@illinois.edu>
#define	MAXCOLOR 128 // available colors: 0 ~ (MAXCOLOR - 1)

void FirstFit(int m, int nnz, int *csrRowPtr, int *csrColInd, Worklist &inwl, int *coloring) {	
	unsigned start = inwl.start;
	unsigned end = inwl.end;
#ifndef ENABLE_OPENMP
	int *forbiddenColors = (int *) malloc(m * sizeof(int));
	for(int i = 0; i < m; i ++) forbiddenColors[i] = m + 1;
#else
	int **forbiddenColors = (int **) malloc(num_omp_threads*sizeof(int*));
	for (int i = 0; i < num_omp_threads; i++) {
		forbiddenColors[i] = (int *) malloc((MAXCOLOR+1)*sizeof(int));
		for(int j = 0; j < MAXCOLOR; j++) forbiddenColors[i][j] = m + 1;
	}
	#pragma omp parallel for
#endif
	for (int i = start; i < end; i++) {
#ifdef ENABLE_OPENMP
		int tid = omp_get_thread_num();
		int vertex = inwl.getItem(i);
#else
		int vertex = i;
#endif
		int row_begin = csrRowPtr[vertex];
		int row_end = csrRowPtr[vertex + 1];
		for (int offset = row_begin; offset < row_end; offset++) {
			int neighbor = csrColInd[offset];
			int color = coloring[neighbor];
#ifdef ENABLE_OPENMP
			forbiddenColors[tid][color] = vertex;
#else
			forbiddenColors[color] = vertex;//forbid this color
#endif
		}
		int vertex_color = 0;
#ifdef ENABLE_OPENMP
		while (vertex_color < MAXCOLOR && forbiddenColors[tid][vertex_color] == vertex)
#else
		while (vertex_color < MAXCOLOR && forbiddenColors[vertex_color] == vertex)
#endif
			vertex_color++;
		assert(vertex_color < MAXCOLOR);
		coloring[vertex] = vertex_color;
	}
}

void conflictDetect(int m, int nnz, int *csrRowPtr, int *csrColInd, Worklist &inwl, Worklist &outwl, int *coloring) {
	unsigned start = inwl.start;
	unsigned end = inwl.end;
#ifdef ENABLE_OPENMP	
	#pragma omp parallel for
#endif
	for (int i = start; i < end; i++) {
		int vertex = inwl.getItem(i);
		int neighbor_offset = csrRowPtr[vertex];
		int num_neighbors = csrRowPtr[vertex + 1] - neighbor_offset;
		for (int j = 0; j < num_neighbors; j++) {
			int neighbor = csrColInd[neighbor_offset + j];
			if (coloring[vertex] == coloring[neighbor] && vertex < neighbor) {
				outwl.push(vertex);
				break;
			}
		}
	}
}

void findMax(int *coloring, int n, int *ncolors) {
	int i;
	for (i = 0; i < n; i++) {
		if (coloring[i] > *ncolors)
			*ncolors = coloring[i];
	}
	*ncolors ++;
}

void color(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring) {
	Worklist inwl, outwl, *inwlptr, *outwlptr, *tmp;
	double starttime, endtime;
	double runtime[ITERATIONS];
	int colors[ITERATIONS];
	int iteration[ITERATIONS];
	for(int i = 0; i < ITERATIONS; i ++) {
		colors[i] = 0;
		iteration[i] = 0;
	}
	inwl.ensureSpace(m);
	outwl.ensureSpace(m);
	for (int i = 0; i < ITERATIONS; i++) {
		inwlptr = &inwl;
		outwlptr = &outwl;
		starttime = rtclock();
		unsigned *range = (unsigned *)malloc(m * sizeof(unsigned));
		for (unsigned j = 0; j < m; j++)
			range[j] = j;
		inwl.pushRange(range, m);
		unsigned wlsz = inwl.getSize();
#ifdef ENABLE_OPENMP
		while (wlsz) {
			++iteration[i];
			//printf("iteration=%d, %d vertices to process\n", iteration, wlsz);
#endif
			FirstFit(m, nnz, csrRowPtr, csrColInd, *inwlptr, coloring);
#ifdef ENABLE_OPENMP
			__syncthreads();
			conflictDetect(m, nnz, csrRowPtr, csrColInd, *inwlptr, *outwlptr, coloring);
			__syncthreads();
			wlsz = outwlptr->getSize();
			tmp = inwlptr; inwlptr = outwlptr; outwlptr = tmp;
			outwlptr->clear();
		}
#endif
		endtime = rtclock();
		findMax(coloring, m, &colors[i]);
		runtime[i] = (1000.0f * (endtime - starttime));
	}
	double total_time = 0.0;
	int total_colors = 0;
	int total_iterations = 0;
	for (int i = 0; i < ITERATIONS; i++) {
		total_time += runtime[i];
		total_colors += colors[i];
		total_iterations += iteration[i];
		printf("[%d %.2f %d] ", colors[i], runtime[i], iteration[i]);
	}   
	double avg_time = (double)total_time / ITERATIONS;
	double avg_colors = (double)total_colors / ITERATIONS;
	double avg_iterations = (double)total_iterations / ITERATIONS;
	printf("\navg_time %f ms avg_colors %.2f avg_iterations %.2f\n", avg_time, avg_colors, avg_iterations);
}
