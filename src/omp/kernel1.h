/*
#include "worklist.h"

#include <vector>
#include <set>

using namespace std;
*/

#define	MAXCOLOR	128

void FirstFit(int m, int nnz, int *csrRowPtr, int *csrColInd, Worklist &inwl, int *coloring)
{	
	unsigned start, end;
	int ii;

	start = inwl.start;
	end = inwl.end;


	#ifdef ENABLE_OPENMP	
	#pragma omp parallel for
	#endif
	for (ii = start; ii < end; ii++) {
		int  j, node, neighbors, neighbor_j;

		node = inwl.getItem(ii);
		int neighboroffset = csrRowPtr[node];
		neighbors = csrRowPtr[node + 1] - neighboroffset;
	
		unsigned v[MAXCOLOR / 32];
		v[0] = 0xfffffffe;
		for (j = 1; j < MAXCOLOR / 32; j++)
			v[j] = 0xffffffff;	

		for (j = 0; j < neighbors; j++) {
			neighbor_j = csrColInd[neighboroffset + j];
			int color_j = coloring[neighbor_j];
			if (color_j)
				v[color_j / 32] &= ~(1 << (color_j % 32));			
		}
		
		int c = 32;
                for (int i = 0; i < MAXCOLOR / 32; i++) {
                        if (v[i] != 0) {
                                v[i] &= -(signed)v[i];
                                if (v[i]) c--;
                                if (v[i] & 0x0000ffff) c -= 16;
                                if (v[i] & 0x00ff00ff) c -= 8;
                                if (v[i] & 0x0f0f0f0f) c -= 4;
                                if (v[i] & 0x33333333) c -= 2;
                                if (v[i] & 0x55555555) c -= 1;
                                break;
                        }
                        else
                                c += 32;
                }
                coloring[node] = c;		
	}
}

void conflictDetect(int m, int nnz, int *csrRowPtr, int *csrColInd, Worklist &inwl, Worklist &outwl, int *coloring)
{
	unsigned start, end;
	int ii;
	//inwl.myItems(start, end);
	start = inwl.start;
	end = inwl.end;
	//printf("inwl=%d, outwl=%d, start=%d, end=%d\n", inwl.getSize(), outwl.getSize(), start, end);

	#ifdef ENABLE_OPENMP	
	#pragma omp parallel for
	#endif
	for (ii = start; ii < end; ii++) {
		int j, node, neighbors, neighbor_j;
		node = inwl.getItem(ii);
		//if (node == -1)
			//continue;
		int neighboroffset = csrRowPtr[node];
		neighbors = csrRowPtr[node + 1] - neighboroffset;
		//neighbors = graph.noutgoing[node];

		for (j = 0; j < neighbors; j++) {
			//neighbor_j = graph.edgessrcdst[graph.psrc[node] + j];
			neighbor_j = csrColInd[neighboroffset + j];
			if (coloring[node] == coloring[neighbor_j] && node < neighbor_j) {
				//printf("c[%d] = c[%d] = %d\n", node, neighbor_j, coloring[node]);
				outwl.push(node);
				break;
			}
		}

		//if (j == neighbors)
			//printf("%d ok\tcolor[%d]=%d\n", node, node, coloring[node]);
	}
}

void findMax(int *coloring, int n, int *ncolors) {
	int i;
	for (i = 0; i < n; i++) {
		//printf("coloring[%d]=%d\n", i, coloring[i]);
		if (coloring[i] > *ncolors)
			*ncolors = coloring[i];
	}
}

void color(int m, int nnz, int *csrRowPtr, int *csrColInd, int *ncolors, int *coloring)
{
	Worklist inwl, outwl, *inwlptr, *outwlptr, *tmp;
	
	double starttime, endtime;
	double runtime;
	
	//int nnodes = graph.nnodes;
	
	inwl.ensureSpace(m);
	outwl.ensureSpace(m);
	inwlptr = &inwl;
	outwlptr = &outwl;
	
	unsigned *range;
	range = (unsigned *)malloc(m * sizeof(unsigned));
	for (unsigned i = 0; i < m; i++)
		range[i] = i;
	//inwl.pushRange(graph.srcsrc, nnodes);
	inwl.pushRange(range, m);

	int iteration = 0;
	unsigned wlsz = inwl.getSize();
	//printf("wlsz=%d, outwl=%d\n", wlsz, outwl.getSize());
	//printf("solving.\n");
		
	starttime = rtclock();	
	#ifdef ENABLE_OPENMP
	while (wlsz) {
		++iteration;
	#endif
	
		//FirstFit(graph, *inwlptr, coloring);
		FirstFit(m, nnz, csrRowPtr, csrColInd, *inwlptr, coloring);
		#ifdef ENABLE_OPENMP
		__syncthreads();
		//printf("ok\n");
		//conflictDetect(graph, *inwlptr, *outwlptr, coloring);
		conflictDetect(m, nnz, csrRowPtr, csrColInd, *inwlptr, *outwlptr, coloring);
		__syncthreads();
		//printf("ok\n");

		//printf("iteration %d:inwl=%d, outwl=%d\n", iteration, wlsz, outwlptr->getSize());
		wlsz = outwlptr->getSize();

		tmp = inwlptr; inwlptr = outwlptr; outwlptr = tmp;
		outwlptr->clear();
	}
		#endif
	endtime = rtclock();
	
	//verify<<<(nnodes - 1) / 1024 + 1, 1024>>>(graph, coloring, correct);
	//CUDA_SAFE_CALL(cudaDeviceSynchronize());
	//if (*correct) {
		//findMax<<<(nnodes - 1) / 1024 + 1, 1024>>>(coloring, nnodes, ncolors);
	findMax(coloring, m, ncolors);
		//CUDA_SAFE_CALL(cudaDeviceSynchronize());
	//}
	
	runtime = (1000.0f * (endtime - starttime));
	printf("runtime=%f\tcolors=%d\t", runtime, *ncolors);
}
