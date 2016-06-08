#ifndef _GRAPHCOLORING_H_ 
#define _GRAPHCOLORING_H_ 

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <math.h> 
#include <cuda_runtime_api.h> 
#include <cuda.h> 
#include <iostream>
using namespace std;


// Should be at least equal to maxDegree of graph + 1
//    if doing that generates an error like: too much local memory, then use commented line 
//    maked OPTION2 instead of OPTION1 in function color & saturation in gaphColoring.cu
//const int TEMP_COLOR_LENGTH = 256;//128; //256;//1024;		
const int TEMP_COLOR_LENGTH = 1000;//128; //256;//1024;		

const int CONFLICT_BLOCK_SIZE = 256;

const int MAXGPUITERATIONS = 50;


#ifdef __cplusplus 
	#define CHECK_EXT extern "C" 
#else 
	#define CHECK_EXT 
#endif 


CHECK_EXT float cudaGraphColoring(int *adjacentList, int *boundaryList, int *graphColors, int *degreeList, 
				int *conflict, int boundarySize, int maxDegree, int graphSize, int & passes, 
				int subsizeBoundary, int _gridSize, int _blockSize, int *startPartitionList, 
				int *endPartitionList, int *randomList, int numRand, int useSDO, int *numOut);


#endif // _GRAPHCOLORING_H_ 

