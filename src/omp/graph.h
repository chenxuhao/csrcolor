#ifndef LSG_GRAPH
#define LSG_GRAPH

#define MYINFINITY	1000000000
#define DISTANCETHRESHOLD	150
#define THRESHOLDDEGREE		10

typedef struct Graph {
	enum {NotAllocated, AllocatedOnHost, AllocatedOnDevice} memory;

	unsigned read(char file[]);
	long unsigned cudaCopy(struct Graph &copygraph);
	unsigned optimize();
	unsigned printStats();
	void     print();

	Graph();
	~Graph();
	unsigned init();
	unsigned allocOnHost();
	unsigned allocOnDevice();
	unsigned dealloc();
	unsigned deallocOnHost();
	unsigned deallocOnDevice();
	unsigned optimizeone();
	unsigned optimizetwo();
	void allocLevels();
	void freeLevels();
	void progressPrint(unsigned maxii, unsigned ii);
	unsigned readFromEdges(char file[]);
	unsigned readFromGR(char file[]);
	unsigned getOutDegree(unsigned src);
	unsigned getDestination(unsigned src, unsigned nthedge);
	unsigned getFirstEdge(unsigned src);
	foru getWeight(unsigned src, unsigned nthedge);

	unsigned nnodes, nedges;
	unsigned *noutgoing, *nincoming, *srcsrc, *psrc, *edgessrcdst;
	foru *edgessrcwt;
	unsigned *levels;
	unsigned source;

	unsigned *maxOutDegree, *maxInDegree;
	unsigned diameter;
	bool foundStats;

} Graph;

unsigned Graph::init() {
	noutgoing = nincoming = srcsrc = psrc = edgessrcdst = NULL;
	edgessrcwt = NULL;
	source = 0;
	nnodes = nedges = 0;
	memory = NotAllocated;

	maxOutDegree = maxInDegree = NULL;
	diameter = 0;
	foundStats = false;

	return 0;
}

unsigned Graph::allocOnHost() {
	edgessrcdst = (unsigned int *)malloc((nedges+1) * sizeof(unsigned int));	// first entry acts as null.
	edgessrcwt = (foru *)malloc((nedges+1) * sizeof(foru));	// first entry acts as null.
	psrc = (unsigned int *)calloc(nnodes+1, sizeof(unsigned int));	// init to null.
	psrc[nnodes] = nedges;	// last entry points to end of edges, to avoid thread divergence in drelax.
	noutgoing = (unsigned int *)calloc(nnodes, sizeof(unsigned int));	// init to 0.
	nincoming = (unsigned int *)calloc(nnodes, sizeof(unsigned int));	// init to 0.
	srcsrc = (unsigned int *)malloc(nnodes * sizeof(unsigned int));

	maxOutDegree = (unsigned *)malloc(sizeof(unsigned));
	maxInDegree = (unsigned *)malloc(sizeof(unsigned));
	*maxOutDegree = 0;
	*maxInDegree = 0;

	memory = AllocatedOnHost;
	return 0;
}

unsigned Graph::deallocOnHost() {
	free(noutgoing);
	free(nincoming);
	free(srcsrc);
	free(psrc);
	free(edgessrcdst);
	free(edgessrcwt);

	free(maxOutDegree);
	free(maxInDegree);
	return 0;
}

unsigned Graph::dealloc() {
	switch (memory) {
		case AllocatedOnHost:
			printf("dealloc on host.\n");
			deallocOnHost();
			break;
		case AllocatedOnDevice:
			printf("dealloc on device.\n");
//			deallocOnDevice();
			break;
	}
	return 0;
}

Graph::Graph() {
	init();
}

Graph::~Graph() {
}

//TODO: make optimizations use the graph api.
unsigned Graph::optimizeone() {
	unsigned int nvv = nnodes;	// no of vertices to be optimized.
	unsigned int insertindex = 1;	// because ii starts with 0.

	for (unsigned ii = 0; ii < nvv; ++ii) {
		unsigned src = srcsrc[ii];
		unsigned dstindex = psrc[src];
		unsigned degree = noutgoing[src];
		if (degree && srcsrc[edgessrcdst[dstindex]] > src + DISTANCETHRESHOLD) {
			unsigned int nee = degree;
			for (unsigned ee = 0; ee < nee; ++ee) {
				unsigned dst = edgessrcdst[dstindex + ee];
				unsigned dstentry = srcsrc[dst];
				// swap insertindex and dst.
				unsigned temp = psrc[insertindex];
				psrc[insertindex] = psrc[dstentry];
				psrc[dstentry] = temp;

				temp = srcsrc[ii];
				srcsrc[ii] = srcsrc[dst];
				srcsrc[dst] = temp;

				if (++insertindex >= nnodes) {
					break;
				}
			}
			if (insertindex >= nnodes) {
				break;
			}
		}
	}
	return 0;
}

unsigned Graph::optimizetwo() {
	// load balance.
	unsigned int nvv = nnodes / 2;
	bool firsthalfsmaller = true;
	unsigned int temp;

	for (unsigned ii = 0; ii < nvv; ++ii) {
		unsigned one = ii;
		unsigned two = nvv + ii;
		unsigned degreeone = noutgoing[one];
		unsigned degreetwo = noutgoing[two];

		if (degreeone > degreetwo && degreeone - degreetwo > THRESHOLDDEGREE && !firsthalfsmaller || degreetwo > degreeone && degreetwo - degreeone > THRESHOLDDEGREE && firsthalfsmaller) {
			temp = srcsrc[one];
			srcsrc[one] = srcsrc[two];
			srcsrc[two] = temp;

			temp = psrc[one];
			psrc[one] = psrc[two];
			psrc[two] = temp;
			firsthalfsmaller = !firsthalfsmaller;
		}
	}
	return 0;
}

unsigned Graph::optimize() {
	optimizeone();
	optimizetwo();
	return 0;
}

void Graph::progressPrint(unsigned maxii, unsigned ii) {
	const unsigned nsteps = 10;
	unsigned ineachstep = (maxii / nsteps);
	if(ineachstep == 0) ineachstep = 1;
	/*if (ii == maxii) {
		printf("\t100%%\n");
	} else*/ if (ii % ineachstep == 0) {
		printf("\t%3d%%\r", ii*100/maxii + 1);
		fflush(stdout);
	}
}

unsigned Graph::readFromEdges(char file[]) {
	std::ifstream cfile;
	cfile.open(file);

	std::string str;
	getline(cfile, str);
	sscanf(str.c_str(), "%d %d", &nnodes, &nedges);

	allocOnHost();
	for (unsigned ii = 0; ii < nnodes; ++ii) {
		srcsrc[ii] = ii;
	}


	unsigned int prevnode = 0;
	unsigned int tempsrcnode;
	unsigned int ncurroutgoing = 0;
	for (unsigned ii = 0; ii < nedges; ++ii) {
		getline(cfile, str);
		sscanf(str.c_str(), "%d %d %d", &tempsrcnode, &edgessrcdst[ii+1], &edgessrcwt[ii+1]);
		if (prevnode == tempsrcnode) {
			if (ii == 0) {
				psrc[tempsrcnode] = ii + 1;
			}
			++ncurroutgoing;
		} else {
			psrc[tempsrcnode] = ii + 1;
			if (ncurroutgoing) {
				noutgoing[prevnode] = ncurroutgoing;
			}
			prevnode = tempsrcnode;
			ncurroutgoing = 1;	// not 0.
		}
		++nincoming[edgessrcdst[ii+1]];

		progressPrint(nedges, ii);
	}
	noutgoing[prevnode] = ncurroutgoing;	// last entries.

	cfile.close();
	return 0;
}

unsigned Graph::readFromGR(char file[]) {
	std::ifstream cfile;
	cfile.open(file);

	// copied from GaloisCpp/trunk/src/FileGraph.h
	int masterFD = open(file, O_RDONLY);
  	if (masterFD == -1) {
	printf("FileGraph::structureFromFile: unable to open %s.\n", file);
	return 1;
  	}

  	struct stat buf;
	int f = fstat(masterFD, &buf);
  	if (f == -1) {
    		printf("FileGraph::structureFromFile: unable to stat %s.\n", file);
    		abort();
  	}
  	size_t masterLength = buf.st_size;

  	int _MAP_BASE = MAP_PRIVATE;
//#ifdef MAP_POPULATE
//  _MAP_BASE  |= MAP_POPULATE;
//#endif

  	void* m = mmap(0, masterLength, PROT_READ, _MAP_BASE, masterFD, 0);
  	if (m == MAP_FAILED) {
    		m = 0;
    		printf("FileGraph::structureFromFile: mmap failed.\n");
    		abort();
  	}

	double starttime, endtime;
	starttime = rtclock();

  	//parse file
  	uint64_t* fptr = (uint64_t*)m;
  	__attribute__((unused)) uint64_t version = le64toh(*fptr++);
  	assert(version == 1);
  	uint64_t sizeEdgeTy = le64toh(*fptr++);
  	uint64_t numNodes = le64toh(*fptr++);
  	uint64_t numEdges = le64toh(*fptr++);
  	uint64_t *outIdx = fptr;
  	fptr += numNodes;
  	uint32_t *fptr32 = (uint32_t*)fptr;
  	uint32_t *outs = fptr32; 
  	fptr32 += numEdges;
  	if (numEdges % 2) fptr32 += 1;
  	unsigned  *edgeData = (unsigned *)fptr32;

	
	// cuda.
	nnodes = numNodes;
	nedges = numEdges;

	printf("nnodes=%d, nedges=%d.\n", nnodes, nedges);
	allocOnHost();

	for (unsigned ii = 0; ii < nnodes; ++ii) {
		// fill unsigned *noutgoing, *nincoming, *srcsrc, *psrc, *edgessrcdst; foru *edgessrcwt;
		srcsrc[ii] = ii;
		if (ii > 0) {
			psrc[ii] = le64toh(outIdx[ii - 1]) + 1;
			noutgoing[ii] = le64toh(outIdx[ii]) - le64toh(outIdx[ii - 1]);
		} else {
			psrc[0] = 1;
			noutgoing[0] = le64toh(outIdx[0]);
		}
		for (unsigned jj = 0; jj < noutgoing[ii]; ++jj) {
			unsigned edgeindex = psrc[ii] + jj;
			unsigned dst = le32toh(outs[edgeindex - 1]);
			if (dst >= nnodes) printf("\tinvalid edge from %d to %d at index %d(%d).\n", ii, dst, jj, edgeindex);
			edgessrcdst[edgeindex] = dst;
			edgessrcwt[edgeindex] = edgeData[edgeindex - 1];

			++nincoming[dst];
			//if (ii == 194 || ii == 352) {
			//	printf("edge %d: %d->%d, wt=%d.\n", edgeindex, ii, dst, edgessrcwt[edgeindex]);
			//}
		}
		progressPrint(nnodes, ii);
	}

	cfile.close();	// probably galois doesn't close its file due to mmap.

	endtime = rtclock();

	printf("read %lld bytes in %0.2f ms (%0.2f MB/s)\n", masterLength, 1000 * (endtime - starttime), (masterLength / 1048576) / (endtime - starttime));

	return 0;
}
unsigned Graph::read(char file[]) {
	if (strstr(file, ".edges")) {
		return readFromEdges(file);
	} else if (strstr(file, ".gr")) {
		return readFromGR(file);
	}
	return 0;
}

unsigned Graph::getOutDegree(unsigned src) {
	if (src < nnodes) {
		return noutgoing[src];
	}
	return 0;
}

unsigned Graph::getDestination(unsigned src, unsigned nthedge) {
	if (src < nnodes && nthedge < getOutDegree(src)) {
		unsigned edge = getFirstEdge(src) + nthedge;
		if (edge && edge < nedges + 1) {
			return edgessrcdst[edge];
		}
		return nnodes;
	}
	if (src < nnodes) {
		printf("Error: %s(%d): node %d: edge %d out of bounds %d.\n", __FILE__, __LINE__, src, nthedge, getOutDegree(src));
	} else {
		printf("Error: %s(%d): node %d out of bounds %d.\n", __FILE__, __LINE__, src, nnodes);
	}
	return nnodes;
}

foru Graph::getWeight(unsigned src, unsigned nthedge) {
	if (src < nnodes && nthedge < getOutDegree(src)) {
		unsigned edge = getFirstEdge(src) + nthedge;
		if (edge && edge < nedges + 1) {
			return edgessrcwt[edge];
		}
		return MYINFINITY;
	}
	if (src < nnodes) {
		printf("Error: %s(%d): node %d: edge %d out of bounds %d.\n", __FILE__, __LINE__, src, nthedge, getOutDegree(src));
	} else {
		printf("Error: %s(%d): node %d out of bounds %d.\n", __FILE__, __LINE__, src, nnodes);
	}
	return MYINFINITY;
}

unsigned Graph::getFirstEdge(unsigned src) {
	if (src < nnodes) {
		unsigned srcnout = getOutDegree(src);
		if (srcnout > 0 && srcsrc[src] < nnodes) {
			return psrc[srcsrc[src]];
		}
		printf("Error: %s(%d): edge %d out of bounds %d.\n", __FILE__, __LINE__, 0, srcnout);
		return 0;
	}
	printf("Error: %s(%d): node %d out of bounds %d.\n", __FILE__, __LINE__, src, nnodes);
	return 0;
}

#endif
