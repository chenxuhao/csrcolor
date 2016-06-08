// Graph coloring 
#include <ctime> 
#include <stdio.h>  
#include <stdlib.h>  
#include <time.h>  
#include <iostream> 
#include <math.h>  
#include <set> 
#include <assert.h>
#include <fstream>

#include <set>

#include "graphColoring.h" 
#include "tree.h"
using namespace std;  


int *degreeList;
int *saturationList;


//----------------------- Utilities -----------------------//
int findPower(int x){
	int num = 2;
	int powerIndex = 1;
	
	while (num <= x){
		powerIndex++;
		num = pow(2,powerIndex);
	}
	
	cout << "Closest power: " << num << endl;
	return num;
}


int findMultiple(int multipleOf, int x){
	int base = multipleOf;
	int powerIndex = 0;
	int num = 0;
	
	while (num <= x){
		powerIndex++;
		num = base * powerIndex;	
	}
	
	//cout << "Closest multiple of " << base << ": " << num << endl;
	return num;
}


int inline min(int n1, int n2) 
{ 
	if (n1>=n2) 
		return n2; 
	else 
		return n1; 
} 




//----------------------- Display -----------------------//
// Author: Pascal
// Displays an adjacencyList
void displayAdjacencyList(int *adjacencyList, int graphSize, int maxDegree){
	//cout << endl << "Adjacency List:" << endl;
	FILE *fp = fopen("adjacencyList.txt", "w");
	for (int i=0; i<graphSize; i++){
		//cout << i << ": ";
		fprintf(fp, "%d: ", i);
		
		for (int j=0; j<maxDegree; j++){
			if (adjacencyList[i*maxDegree + j] != -1)
				//cout << adjacencyList[i*maxDegree + j] << " ";
				fprintf(fp, "%d ", adjacencyList[i*maxDegree + j]);
			else 
				break;
		}
		
		//cout << endl;
		fprintf(fp, "\n");
	}
}





//----------------------- Graph initializations -----------------------//
// Author: Pascal & Shusen
// reads a matrix market format into an adjacency matrix
void readGraph(int *&adjacencyMatrix, const char *filename, int _gridSize, int _blockSize, int &graphSizeRead, int &graphSize, long &edgeSize){
	char comments[512], weightedAns;
	int graphSizeX, graphSizeY, from, to, numEdges, weightedGraph;
	float weight;
	
	
	
	numEdges = 0;
	
	ifstream graphFile(filename);
	
	
	if (graphFile.is_open())
	{
		while (graphFile.peek() == '%'){
			graphFile.getline(comments,512);
		}
		
		graphFile >> graphSizeX >> graphSizeY >> edgeSize;
		cout << "Rows: " << graphSizeX << "  ,  Columns: " <<  graphSizeY << "   - Number of edges: " << edgeSize << endl;
		
		if (graphSizeX != graphSizeY){
			cout << "Non Symmetric graph!" << endl;
			exit(1);
		}
		else 
		{
			cout << "Is it a weighted graph(y/n): ";
			cin >> weightedAns;

			if (weightedAns == 'y')
				weightedGraph = 1;
			else
				weightedGraph = 0;
			
			
			graphSizeRead = graphSizeX;
			graphSize = findMultiple(_gridSize*_blockSize, graphSizeRead);
			
			adjacencyMatrix = new int[graphSize * graphSize];
			memset(adjacencyMatrix, 0, graphSize * graphSize *sizeof(int));
			
			for (int i=0; i<edgeSize; i++){
				if (weightedGraph == 1)
					graphFile >> from >> to >> weight;	
				else
					graphFile >> from >> to;
				
				if (!(from == to)){
					numEdges++;
					adjacencyMatrix[(from-1)*graphSize + (to-1)] = 1;
					adjacencyMatrix[(to-1)*graphSize + (from-1)] = 1;
					
					/*
					 if (weightedGraph == 1)
					 cout << from << " , " << to << " : " << weight << endl; 
					 else
					 cout << from << " , " << to << endl;
					 */
				}
			}
		}
	}
	else {
		cout << "Reading " << filename << " failed!" << endl;
		exit(1);
	}
	
	edgeSize = numEdges;
	cout << "Graph: " << graphSizeRead << " - " <<  graphSize << " - " <<  edgeSize << endl;
	cout << "File " << filename << " was successfully read!" << endl;
}



// Author: Pascal & Shusen
// reads a matrix market format into an adjacency list
void readGraphAdj(int *&adjacencyList, const char *filename, int _gridSize, int _blockSize, int &graphSizeRead, int &graphSize, long &edgeSize, int &maxDegree, char weighted, int interactive, bool useMetis)
{		
	char comments[512], weightedAns;
	int graphSizeX, graphSizeY, numUsefulEdges, from, to, weightedGraph;
	int *edgesCount, *degreeList;
	double weight;
	
	numUsefulEdges = maxDegree = 0;
	
	
	ifstream graphFile(filename);
	
	
	if (graphFile.is_open())
	{
		//
		// Pass 1: get the degree of each node
		//
		while (graphFile.peek() == '%'){
			graphFile.getline(comments,512);
		}
		
		graphFile >> graphSizeX >> graphSizeY >> edgeSize;
		if (interactive != 2)
			cout << "Rows: " << graphSizeX << "  ,  Columns: " <<  graphSizeY << "   - Number of edges: " << edgeSize << endl;
		
		if (graphSizeX != graphSizeY){
			cout << "Non Symmetric graph!" << endl;
			exit(1);
		}
		else 
		{
			if (weighted == 'i'){
				cout << "Is it a weighted graph(y/n): ";
			//if (weighted == 'i')
				cin >> weightedAns;
			}
			else
				weightedAns = weighted;

			
			if (weightedAns == 'y')
				weightedGraph = 1;
			else
				weightedGraph = 0;
			
			graphSizeRead = graphSizeX;
			if (useMetis == false)
				graphSize = findMultiple(_gridSize*_blockSize, graphSizeRead);
			else
				graphSize = graphSizeX;		
	
			degreeList = new int[graphSize];
			memset(degreeList, 0, graphSize*sizeof(int));
			
			
			for (int i=0; i<edgeSize; i++){
				if (weightedGraph == 1){
					graphFile >> from >> to >> weight;	
					//	cout << from << " , " << to << " : " << weight << endl;
				}
				else{
					graphFile >> from >> to;	
					//	cout << from << " , " << to << endl;
				}
				
				
				if (from != to){
					degreeList[(from-1)]+= 1;
					degreeList[(to-1)]+= 1;
					
					numUsefulEdges++;	/* 可能重复计算 */
				}
			}
		}
		printf("pass 1 ok\n");		
		
		//
		// Get max degree and reset
		//
		for (int i=0; i<graphSize; i++){
			//cout << i << ": " << degreeList[i] << endl;  
			if (maxDegree < degreeList[i]){
				maxDegree = degreeList[i];
				//				cout << i << "  Max degree: " << maxDegree << endl;  
			}
		}
		
		
		edgesCount = new int[graphSize];
		memset(edgesCount, 0, graphSize*sizeof(int));
		
		adjacencyList = new int[graphSize*maxDegree];
		memset(adjacencyList, -1, graphSize*maxDegree*sizeof(int));
		
		printf("get adjList ok\n");
		
		graphFile.seekg(0, ios::beg);
		
		
		
		//
		// Pass 2: Fill in the adjacency List
		//
		while (graphFile.peek() == '%'){
			graphFile.getline(comments,512);
		}
		
		
		graphFile >> graphSizeX >> graphSizeY >> edgeSize;
		//cout << graphSizeX << " " <<  graphSizeY << " - " << numEdges << endl;
		
		if (graphSizeX != graphSizeY){
			cout << "Non Symmetric graph!" << endl;
			exit(1);
		}
		else 
		{
			for (int i=0; i<edgeSize; i++){
				
				if (weightedGraph == 1){
					graphFile >> from >> to >> weight;	
					//cout << from << " , " << to << " : " << weight << endl;
				}
				else{
					graphFile >> from >> to;	
					//cout << from << " , " << to << endl;
				}
				
				
				
				if (from != to){
					adjacencyList[(from-1)*maxDegree + edgesCount[(from-1)]] = (to-1);
					edgesCount[(from-1)]++;
					
					
					adjacencyList[(to-1)*maxDegree + edgesCount[(to-1)]] = (from-1);
					edgesCount[(to-1)]++;
				}
			}
		}
		printf("pass 2 ok\n");
	}
	else {
		cout << "Reading " << filename << " failed!" << endl;
		exit(1);
	}
	
	graphFile.close();
	
	
	edgeSize = numUsefulEdges;
	if (interactive != 2){
		cout << "Graph: " << graphSizeRead << " - " <<  graphSize << " - " <<  edgeSize << endl;
		cout << "Max degree: " << maxDegree << endl;
		cout << "File " << filename << " was successfully read!" << endl;
	}

	delete []degreeList;
	delete []edgesCount;
}



// Author: Pascal 
// Genetates a graph 
void generateMatrix(int *matrix, int graphSize, int num){  
	int x, y;  
	int count = 0;  
	
	while (count < num){  
		x = rand()%graphSize;  
		y = rand()%graphSize;  
		
		if (matrix[x*graphSize + y] != 1)          // if not already assigned an edge 
			if (x != y){  
				matrix[x*graphSize + y] = 1;       // non directional graph  
				matrix[y*graphSize + x] = 1;  
				count++;  
			}  
	}  
}  



// Author:Peihong
// node index start from 1
// gets an adjacency list from an adjacencyMatrix
void getAdjacentList(int *adjacencyMatrix, int *adjacentList, int size, int maxDegree)
{
	for (int i=0; i<size; i++){ 
		int nbCount = 0;
		for (int j=0; j<size; j++){
			if ( adjacencyMatrix[i*size + j] == 1)  
			{
				adjacentList[i*maxDegree + nbCount] = j;
				nbCount++;
			}
		}
	}
}



// Author: Pascal 
// get the degree information for a graph 
int getMaxDegree(int *adjacencyMatrix, int size, int &avgDegree){  
	int maxDegree = 0;   
	int degree;  
	avgDegree = 0;
	
	for (int i=0; i<size; i++){  
		degree = 0;  
		
		for (int j=0; j<size; j++){         
			if (adjacencyMatrix[i*size + j] == 1){
				degree++; 
				//	cout << i << " , " << j << endl;
			}
			
		}
		
		if (degree > maxDegree)  
			maxDegree = degree;  
		
		avgDegree += degree;
	}  
	
	avgDegree = avgDegree/size;
	
	return maxDegree;  
}  



// Author: Pascal
// get the degree of each element in the graph and returns the maximum degree
int getDegreeList(int *adjacencyList, int *degreeList, int sizeGraph, int maxDegree){
	int avgDegree = 0;
	
    for (int i=0; i<sizeGraph; i++){
        int count = 0;
        
        for (int j=0; j<maxDegree; j++){
            if (adjacencyList[i*maxDegree + j] != -1)
                count++;
            else
                break;  
        }
		
        degreeList[i] = count;
		avgDegree += count;
    }
	
	return (avgDegree/sizeGraph);
}



// Author: Peihong & Pascal
int getBoundaryList(int *adjacencyList, int *boundaryList, int graphSize, int maxDegree, int _gridSize, int _blockSize, int *startPartitionList, int *endPartitionList){  
	int boundaryCount = 0;
	set<int> boundarySet; 
	boundarySet.clear(); 
	
	int start, end, node;
	int partitionIndex = 0; 
	//int subSize = graphSize/(_gridSize*_blockSize);
	//cout << "SubSize = " << subSize << endl;
	
	for (int i=0; i<graphSize; i++)
	{  
		//int subIdx = i/(float)subSize;
		//start = subIdx * subSize;
		//end = min( (subIdx + 1)*subSize, subSize );
		
		if (!(i < endPartitionList[partitionIndex]))
			partitionIndex++;
		
		start = startPartitionList[partitionIndex];
		end = endPartitionList[partitionIndex];
		
		
		
		for (int j=0; j<maxDegree; j++){ 
			node = adjacencyList[i*maxDegree + j]; 
			if (node == -1) 
				break;
			else
				if ((node < start) || (node >= end)){
					boundarySet.insert(i);
					break;		//	we just need one proof of that
				}
		}
		
	} 
	
	boundaryCount = boundarySet.size();
	
	set<int>::iterator it = boundarySet.begin(); 
	for (int i=0; it != boundarySet.end(); it++)  
	{ 
		boundaryList[i] = *it;
		i++; 
	}  
	
	return boundaryCount;  
} 

/* 8个参赛的getBoundaryNodes */
int getBoundaryNodes(int *adjacencyList, int *boundaryList, int graphSize, int maxDegree, int _gridSize, int _blockSize, int *startPartitionList, int *endPartitionList){
	int boundaryCount = 0;
	int boundaryIndex = 0;
	int partitionIndex = 0;
	int start, end, node;
	
	int *tempBoundaries;
	tempBoundaries = new int[graphSize];	/* 存储边界结点 */
	memset(tempBoundaries,0,graphSize*sizeof(int));
	
	for (int i=0; i<graphSize; i++){
		//while (endPartitionList[partitionIndex] < i)
		//	partitionIndex++;
		
		if (!(i < endPartitionList[partitionIndex]))
			partitionIndex++;	
		
		start = startPartitionList[partitionIndex];
		end = endPartitionList[partitionIndex];
		
		for (int j=0; j<maxDegree; j++){
			
			node = adjacencyList[i*maxDegree + j]; 
			if (node == -1)
				break;
			else
				if ((node < start) || (node >= end)){
					tempBoundaries[i] = 1;
					boundaryCount++;
					break;
				} 
		}
	}
	
	
	for (int i=0; i<graphSize; i++){
		if (tempBoundaries[i] == 1){
			boundaryList[boundaryIndex] = i;
			boundaryIndex++;
		}
		
		/* debug */
		//else
			//printf("node %d is not boundary node\n", i);
	}
	
	/* debug */
	printf("boundaryCount=%d\tgraphSize=%d\n", boundaryCount, graphSize);
	
	
	delete []tempBoundaries;
	return boundaryCount;
}

/* 9个参赛的getBoundaryNodes */
int getBoundaryNodes(int *adjacencyList, int *boundaryList, int *numOutside, int graphSize, int maxDegree, int _gridSize, int _blockSize, int *startPartitionList, int *endPartitionList){
	int boundaryCount = 0;
	int boundaryIndex = 0;
	int partitionIndex = 0;
	int start, end, node;
	int hasOutside = 0;
	
	int *tempBoundaries;
	tempBoundaries = new int[graphSize];
	memset(tempBoundaries,0,graphSize*sizeof(int));
	
	for (int i=0; i<graphSize; i++){
		//while (endPartitionList[partitionIndex] < i)
		//      partitionIndex++;
		
		if (!(i < endPartitionList[partitionIndex]))
			partitionIndex++;
		
		start = startPartitionList[partitionIndex];
		end = endPartitionList[partitionIndex];
		
		hasOutside = 0;
		for (int j=0; j<maxDegree; j++){
			
			node = adjacencyList[i*maxDegree + j];
			if (node == -1)
				break;
			else
				if ((node < start) || (node >= end)){
					tempBoundaries[i] = 1;
					if (hasOutside == 0){
						boundaryCount++;
						hasOutside = 1;
						numOutside[i]++;
					}else
						numOutside[i]++;
				}
		}
	}
	
	
	for (int i=0; i<graphSize; i++){
		if (tempBoundaries[i] == 1){
			boundaryList[boundaryIndex] = i;
			boundaryIndex++;
		}
	}
	
	
	delete []tempBoundaries;
	return boundaryCount;
}



//----------------------- Fast Fit Graph Coloring -----------------------//
// Author: Pascal & Shusen
// GraphColor Adjacency list
int colorGraph_FF(int *adjList, int *colors, int size, int maxDegree){  
	int numColors = 0;  
	int i, j;  
	
	int *tempColors;  
	tempColors = new int[maxDegree+1];  
	
	
	for (i=0; i<size; i++)  
	{                 
		// initialize degree array  
		memset(tempColors, 0, (maxDegree+1)*sizeof(int)); 
		
		
		// check the colors  
		for (j=0; j<maxDegree; j++)
			if (adjList[i*maxDegree + j] == -1)
				break;
			else
				tempColors[ colors[adjList[i*maxDegree + j]] ] = 1;   // set connected spots to 1
		
		
		for (j=1; j<=maxDegree; j++)  
			if (tempColors[j] != 1){  
				colors[i] = j;  
				break;  
			}  
		
		if (colors[i] > numColors)  
			numColors=colors[i];  
	}  
	
	delete[] tempColors; 
	
	return numColors;  
}  




//----------------------- SDO Improved Graph Coloring -----------------------//
// Author: Pascal
// returns the degree of that node
int degree1(int vertex, int *degreeList){
	return degreeList[vertex];
}



// Author: Pascal
// return the saturation of that node
int saturation(int vertex, int *adjacencyList, int *graphColors, int maxDegree){
    int saturation = 0;
    int *colors = new int[maxDegree+1];
	
    memset(colors, 0, (maxDegree+1)*sizeof(int));           // initialize array
	
	
    for (int i=0; i<maxDegree; i++){
        if (adjacencyList[vertex*maxDegree + i] != -1)
			//  colors[ graphColors[vertex] ] = 1;                      // at each colored set the array to 1
			colors[ graphColors[adjacencyList[vertex*maxDegree + i]] ] = 1;                      // at each colored set the array to 1
        else
            break;
    }
	
	
    for (int i=1; i<maxDegree+1; i++)                                       // count the number of 1's but skip uncolored
        if (colors[i] == 1)
            saturation++;
	
	delete[] colors; 	
	
    return saturation;
}




// Author: Pascal
// colors the vertex with the min possible color
int color(int vertex, int *adjacencyList, int *graphColors, int maxDegree, int numColored){
    int *colors = new int[maxDegree + 1];
    memset(colors, 0, (maxDegree+1)*sizeof(int));   
    
    if (graphColors[vertex] == 0)
        numColored++;
	else{
		//cout << "Recoloring" << endl;
		graphColors[vertex] = 0;
	}
	//	else
	//		cout << "Old color: " << graphColors[vertex] << "   ";
    
    for (int i=0; i<maxDegree; i++)                                         // set the index of the color to 1
        if (adjacencyList[vertex*maxDegree + i] != -1)
            colors[  graphColors[  adjacencyList[vertex*maxDegree + i]  ]  ] = 1;
        else {
            break;
        }
	
    
	
    for (int i=1; i<maxDegree+1; i++)                                       // nodes still equal to 0 are unassigned
        if (colors[i] != 1){
            graphColors[vertex] = i;
	    //cout << vertex << " : " << i << endl;
            break;
        }
    
	delete[] colors; 
	
    return numColored;
}



int sdoImO(int *adjacencyList, int *graphColors, int *degreeList, int sizeGraph, int maxDegree){
    int satDegree, numColored, max, index;
    numColored = 0;
    int iterations = 0;
    
	cout << endl << "SDO Improved - Allocation outside"; 
	
    while (numColored < sizeGraph){
        max = -1;
        
        for (int i=0; i<sizeGraph; i++){
            if (graphColors[i] == 0)                        // not colored
            {
                satDegree = saturation(i,adjacencyList,graphColors, maxDegree);
				
                if (satDegree > max){
                    max = satDegree;
                    index = i;                              
                }
				
                if (satDegree == max){
                    if (degree1(i,degreeList) > degree1(index,degreeList))
                        index = i;
                }
			}
		}

//		cout << "Num colored: " << numColored << endl;				
		numColored = color(index,adjacencyList,graphColors, maxDegree, numColored);
		iterations++;
    }
    
    return iterations;
}



int sdoIm(int *adjacencyList, int *graphColors, int *degreeList, int sizeGraph, int maxDegree){
    int satDegree, numColored, max, index;
    numColored = 0;
    int iterations = 0;

    cout << endl << "Normal SDO Improved!";

    while (numColored < sizeGraph){
        max = -1;

        for (int i=0; i<sizeGraph; i++){
            if (graphColors[i] == 0)                        // not colored
            {
                satDegree = saturation(i,adjacencyList,graphColors, maxDegree);

                if (satDegree > max){
                    max = satDegree;
                    index = i;
                }

                if (satDegree == max){
                    if (degree1(i,degreeList) > degree1(index,degreeList))
                        index = i;
                }
			}

//          cout << "Num colored: " << numColored << endl;                          
            //if (graphColors[index] == 0){
            	numColored = color(index,adjacencyList,graphColors, maxDegree, numColored);
                iterations++;
            //}
		}
    }

    return iterations;
}


// Author: Pascal
// main driver function for graph coloring
int sdoIm(int *adjacencyList, int *graphColors, int *degreeList, int sizeGraph, int maxDegree, int numPartitions, int *startPartitionList, int *endPartitionList ){
    int satDegree, numColored, max, index;
    numColored = 0;
    int iterations = 0;
    int start, end, subGraphSize;

	cout << endl << "Parts SDO - poor color quality!";

	for (int k=0; k<numPartitions; k++){
		
		start = startPartitionList[k];
		end = endPartitionList[k];
		subGraphSize = end - start;
 
		numColored = 0;
    	while (numColored < subGraphSize){
        	max = -1;
        
        	for (int i=start; i<end; i++){
			//for (int i=0; i<sizeGraph; i++){
            	if (graphColors[i] == 0)                        // not colored
            	{
                	satDegree = saturation(i,adjacencyList,graphColors, maxDegree);
				
                	if (satDegree > max){
                    	max = satDegree;
                    	index = i;                              
                	}
				
                	if (satDegree == max){
                    	if (degree1(i,degreeList) > degree1(index,degreeList))
                        	index = i;
                	}	
				
					numColored = color(index,adjacencyList,graphColors, maxDegree, numColored);
					iterations++;
           	 	}	 
			
				//numColored = color(index,adjacencyList,graphColors, maxDegree, numColored);
				//iterations++;
        	}
    	}
	}

    return iterations;
}





//----------------------- New CPU representation -----------------------//





struct classcomp {
	bool operator() (const int& lhs, const int& rhs) const
	{
		if (saturationList[rhs] > saturationList[lhs])
			return true;
		else 
			if (saturationList[rhs] == saturationList[lhs])
				if (degreeList[rhs] > degreeList[lhs])
					return true;
				else 
					if (degreeList[rhs] == degreeList[lhs])
						if (rhs > lhs)
							return true;
		
		return false;
	}
	
};

void color(int vertex, int *adjacencyList, int *graphColors, int *colorList, int maxDegree, multiset<int, classcomp> &treeGraph){
	int colorAssigned = -1;
	int neighbour;
	int neighTemp;
	
	multiset<int, classcomp>::iterator it;
	
	//Color the node
	for (int i=1; i<maxDegree; i++)
		if (colorList[vertex*maxDegree + i] == 0){
			colorAssigned = i;
			graphColors[vertex] = colorAssigned;
			break;
		}
	
	//cout << vertex << " : " << saturationList[vertex] << " , " <<degreeList[vertex] << " - " << colorAssigned << endl;
	
	for (int i=0; i<degreeList[vertex]; i++){
		neighbour = adjacencyList[vertex*maxDegree + i];
		
		if (graphColors[neighbour] != 0)
			continue;
		
		if (colorList[neighbour*maxDegree + colorAssigned] == 0){
			//remove from graph
			it = treeGraph.find(neighbour);	
			neighTemp = *it;
			treeGraph.erase(it);
			
			// update saturation 
			colorList[neighbour*maxDegree + colorAssigned] = colorAssigned;
			saturationList[neighbour]++;
			
			//reinsert in graph
			treeGraph.insert(neighTemp);
		}
	}
	
}




void seqGraphColoring(int *adjacencyList, int *graphColors, int *degreeList, int maxDegree, int graphSize){
	int index, saturation, degree;
	
	multiset<int, classcomp> treeGraph;
	multiset<int, classcomp>::iterator it;
	
	
	int *colorList = new int[graphSize*maxDegree];  
	memset(colorList, 0, graphSize*maxDegree*sizeof(int));
	
	saturationList = new int[graphSize];  
	memset(saturationList, 0, graphSize*sizeof(int));
	

	
	cout << "Tree SDO New" << endl;
	
	// represent the graph information as a tree sorted by saturation and degree
	for (int i=0; i<graphSize; i++){
		treeGraph.insert(i);
	}
	

	// Color the graph
	for (int i=0; i<graphSize; i++){
		it = treeGraph.end();	it--;
		
		index = *it;
		treeGraph.erase(it);
		
		color(index, adjacencyList, graphColors, colorList, maxDegree, treeGraph);
	}
	
	delete []colorList;
	delete []saturationList;
}






// Author: Pascal
// colors the vertex with the min possible color
/*
void color(int vertex, int *adjacencyList, int *graphColors, int *colorList, int *saturationList, int *degreeList, int maxDegree, tree &graph){
	int colorAssigned = -1;
	int neighbour;
	node* temp;
	
	//Color the node
	for (int i=1; i<maxDegree; i++)
		if (colorList[vertex*maxDegree + i] == 0){
			colorAssigned = i;
			graphColors[vertex] = colorAssigned;
			break;
		}
	
	for (int i=0; i<degreeList[vertex]; i++){
		neighbour = adjacencyList[vertex*maxDegree + i];
		if (graphColors[neighbour] != 0)
			continue;
		
		if (colorList[neighbour*maxDegree + colorAssigned] == 0){
			//remove from graph
			temp = graph.remove(neighbour,saturationList[neighbour],degreeList[neighbour]);
			
			colorList[neighbour*maxDegree + colorAssigned] = colorAssigned;
			saturationList[neighbour]++;
			
			//reinsert in graph
			temp->setSaturation(saturationList[neighbour]);
			temp->setLeft(NULL);
			temp->setRight(NULL);
			graph.insert(temp);
			
		//	graph.displayTreeRML(graph.getTop());
		}
	}
	
}
*/








/*
void seqGraphColoring(int *adjacencyList, int *graphColors, int *degreeList, int maxDegree, int graphSize){
	int index, saturation, degree;
	tree graph;
	
	int *colorList = new int[graphSize*maxDegree];  
	memset(colorList, 0, graphSize*maxDegree*sizeof(int));
	
	int *saturationList = new int[graphSize];  
	memset(saturationList, 0, graphSize*sizeof(int));
	
	
		
	node *nodes = new node[graphSize];

	cout << endl << "Tree SDO";
	
	// represent the graph information as a tree sorted by saturation and degree
	for (int i=0; i<graphSize; i++){
	//	cout << i << " : " << 0 << " , " <<  degreeList[i] << endl;
		nodes[i].setKSD(i, 0, degreeList[i]);
		graph.insert(&nodes[i]);
	}
	
	//cout << "Graph allocation successful" << endl;
	//graph.displayTreeRML(graph.getTop());
	
	
	// Color the graph
	for (int i=0; i<graphSize; i++){
	//	cout << endl << "Graph: " << endl;	graph.displayTreeRML(graph.getTop());
		
		graph.findBiggest(index, saturation, degree);	// get the node with the highest saturation then degree		
		graph.remove(index, saturation, degree);
		
		color(index, adjacencyList, graphColors, colorList, saturationList, degreeList, maxDegree, graph);
	}
	
	delete []colorList;
	delete []saturationList;
}
*/



//----------------------- Conflict Solve -----------------------//
// Author: Pascal
void conflictSolveSDO(int *adjacencyList, int *conflict, int conflictSize, int *graphColors, int *degreeList, int sizeGraph, int maxDegree){
    int satDegree, numColored, max, index, vertex;
    numColored = 0;
    
	
    while (numColored < conflictSize){
        max = -1;
        
        for (int i=0; i<conflictSize; i++){
			vertex = conflict[i];
			
            if (graphColors[vertex] == 0)                        // not colored
            {
                satDegree = saturation(vertex, adjacencyList, graphColors, maxDegree);
				
                if (satDegree > max){
                    max = satDegree;
                    index = vertex;                              
                }
				
                if (satDegree == max){
                    if (degree1(vertex,degreeList) > degree1(index,degreeList))
                        index = vertex;
                }
            }
			
			if (graphColors[index] == 0)
				numColored = color(index,adjacencyList,graphColors, maxDegree, numColored);
        }
		//numColored = color(index,adjacencyList,graphColors, maxDegree, numColored);
    }
}







// Author: Pascal & Shusen
// Solves conflicts using Fast Fit
/* 串行FF冲突解决 */
void conflictSolveFF(int *Adjlist, int size, int *conflict, int conflictSize, int *graphColors, int maxDegree){
	int i, j, vertex, *colorList;
	colorList = new int[(maxDegree+1)];	/* 保存所有邻接点的颜色 */
	
	/* 遍历所有冲突结点 */
	for (i=0; i<conflictSize; i++){
		memset(colorList, 0, (maxDegree+1)*sizeof(int));
		
        vertex = conflict[i];
        
		/* 保存所有邻接点的颜色 */
        for (j=0; j<maxDegree; j++){                                            // cycle through the graph
			if ( Adjlist[vertex*maxDegree + j] != -1 )                      		//      check if node is connected
				colorList[ graphColors[Adjlist[vertex*maxDegree + j]] ] = 1;
			else 
                break;    
        }
		
		/* 给冲突结点着最小的可以着的颜色 */
        for (j=1; j<=maxDegree; j++){                                       	// check the colorList array
			if (colorList[j] != 1){                                         //    at the first spot where we have a color not assigned
				graphColors[vertex] = j;                         //       we assign that color to the node and
                break;                                                      //   	 exit to the next
            }
        }
		
	}
	
	delete []colorList; 
}




//----------------------- Metis -----------------------//
// Metis related stuff...
// Author: Pascal & Shusen
// Creates outpuf for metis file
void createMetisInput(int *adjacencyList, int graphSize, int numEdges, int maxDegree, string metisInputFile, int interactive){
	string metisInputFilename;

	if (interactive != 2){	
		cout << endl << "Creating file to send to metis ..." << endl;
		cout << "Enter filename for file: ";
	}
	if (interactive == 1)
		cin >> metisInputFilename;
	else{
		metisInputFilename = metisInputFile;
		//cout << metisInputFilename << endl;
	}
	
	if (interactive != 2){	
		cout << "Graph Size: " << graphSize << endl;
		cout << "Num Edges: " << numEdges << endl;
	}
	
	ofstream myfile (metisInputFilename.c_str());
  	if (myfile.is_open())
  	{
		myfile << graphSize << " " << numEdges << "\n";
		
		for (int i=0; i<graphSize; i++){ 
			
			for (int j=0; j<maxDegree; j++)  
				if (adjacencyList[i*maxDegree + j] == -1)
					break;
				else
					myfile << (adjacencyList[i*maxDegree + j]) + 1 << " ";
			
			myfile << endl;  
		}  
		
		myfile.close();
  	}
  	else {
		cout << "Unable to open file to write";
		exit(0);
	}
}



// Author: Pascal & Peihong 
// Reads in metis partitioned file
void readMetisOutput(int *partitionList, int graphSize, string metisOutputFile, int interactive){
	string metisOutputFilename;
	
	if (interactive != 2){
		cout << endl << "Reading partitioned metis file..." << endl;
		cout << "Enter filename to read from (e.g. metisInput2048.txt.part.256): ";
	}
	if (interactive == 1)
		cin >> metisOutputFilename;
	else{
		metisOutputFilename = metisOutputFile;
		if (interactive != 2)
			cout << metisOutputFilename << endl;
	}
	
	
	ifstream metisFile(metisOutputFilename.c_str());
	if (metisFile.is_open()){
		for(int i=0; i<graphSize; i++)
		{
			metisFile >> partitionList[i];
		}
		metisFile.close();
		
		/*
		 cout << "Partition List" << endl;
		 for(int i=0; i<graphSize; i++)
		 cout << i << ":: " << partitionList[i] << endl;
		 */
		
	}
	else {
		cout << "Reading in file failed" << endl;
		exit(0);
	}
}



// Author: Pascal & Pihong
// Output file for use by metis
// Input partitioned file
// gets the new adjacency list
int metis(int *adjacencyList, int *newAdjacencyList, int graphSize, int numEdges, int maxDegree, int *startPartitionList, int *endPartitionList, int numMetisPartitions, string metisInputFile, string metisOutputFile, int interactive){
	int *partitionList = new int[graphSize];
	int *newGraphOrdering = new int[graphSize];
	
	int *adjacencyListOrg = new int[graphSize*maxDegree];	// created so as not to modify the original adjacency List
	memcpy(adjacencyListOrg, adjacencyList, graphSize*maxDegree*sizeof(int));
	
	
	int currentPosition = 0;
	int numPartitions = 256;
	
	memset(newAdjacencyList, -1, graphSize*maxDegree*sizeof(int));
	
	
	
	createMetisInput(adjacencyListOrg, graphSize, numEdges, maxDegree, metisInputFile, interactive);

	if (interactive != 2)	
		cout << endl << "Enter the number of partitions used in metis: ";
	
	if (interactive == 1)
		cin >> numPartitions;
	else
		numPartitions = numMetisPartitions;
	
	
	readMetisOutput(partitionList, graphSize, metisOutputFile, interactive);
	
	
	
	// Get the maximum and minimum in each partition
	int min, max, count, partitionMin, partitionMax, startPartitionCount ,endPartitionCount;
	min = 1000000;
	max = -1;
	partitionMin = partitionMax = -1;
	startPartitionCount = endPartitionCount = 0;
	
	
	for (int i=0; i<numPartitions; i++){
		//cout << " i " << partitionList[i] << endl;
		count = 0;
		startPartitionCount = endPartitionCount;
		
		for (int j=0; j<graphSize; j++){
			if (partitionList[j] == i){
				count++;
			}

		}
		endPartitionCount += count;
		
		startPartitionList[i] = startPartitionCount;
		endPartitionList[i] = endPartitionCount;
		
		if (count > max){
			max = count;
			partitionMax = i;
		}
		
		if (count < min){
			min = count;
			partitionMin = i;
		}
	
//		cout << i << " : " << partitionList[j] << " - " << count << endl;
	}

	if (interactive != 2){	
		cout << "Min in partiton: " << min << "  for partition: " << partitionMin << endl;
		cout << "Max in parition: " << max << "  for partition: " << partitionMax << endl;
	}
	else{
		cout << min <<  "  "  << max << " "; 
	}
	/*
	 cout << "Partitions list:" << endl;
	 for (int i=0; i<numPartitions; i++)
	 cout << i << "-   start: " << startPartitionList[i] << "    end: " << endPartitionList[i] 
	 <<  "  size: " << (endPartitionList[i] - startPartitionList[i]) << endl;
	 cout << endl;
	 */
	
	
	// Gets the new Ordering of the nodes
	for (int i=0; i<numPartitions; i++)
		for (int j=0; j<graphSize; j++){
			if (partitionList[j] == i){
				newGraphOrdering[j] = currentPosition;
				currentPosition++;
			}
		}
	
	/*
	 cout << "New Ordering" << endl;
	 for (int i=0; i<graphSize; i++)
	 {
	 cout << i << ": " << newGraphOrdering[i] << endl;
	 }
	 */
	
	
	
	// Replaces the nodes in the adjacency list by the new ordering numbers
	for (int i=0; i<graphSize; i++){
		for (int j=0; j<maxDegree; j++){
			int node = adjacencyListOrg[i*maxDegree + j];
			
			if (adjacencyListOrg[i*maxDegree + j] != -1){
				adjacencyListOrg[i*maxDegree + j] = newGraphOrdering[node];
			}
			else
				break;
		}	
	}
	
	
	
	// Places the different lists of the adjacency list in the right place
	for (int i=0; i<graphSize; i++){
		int newPosn = newGraphOrdering[i];
		
		for (int j=0; j<maxDegree; j++){
			if 	(adjacencyListOrg[i*maxDegree + j] != -1)
				newAdjacencyList[newPosn*maxDegree + j] = adjacencyListOrg[i*maxDegree + j];
			else
				break;
		}	
	}
	
	/*
	 cout << "New Adjacency List:" << endl; 
	 displayAdjacencyList(newAdjacencyList, graphSize, maxDegree);
	 */
	
	delete []adjacencyListOrg;
	delete []partitionList;
	delete []newGraphOrdering;
	
	return numPartitions;
}




//----------------------- Checking for error -----------------------//
// Checking if coloring has been done properly...
// Author: Pascal 
// Checks if a graph has conflicts or not from Adjacency Matrix
void checkCorrectColoring(int *adjacencyMatrix, int *graphColors, int graphSize){ 
	int numErrors = 0;
	int maxColor = -1; 
	
	cout << endl << "==================" << endl << "Error checking for Graph" << endl; 
	
	for (int i=0; i<graphSize; i++)                 // we check each row 
	{ 
		int nodeColor = graphColors[i]; 
		int numErrorsOnRow = 0; 
		

		for (int j=0; j<graphSize;j++){ // check each column in the matrix 
			
			// skip itself 
			if (i == j) 
				continue; 
			
			if (adjacencyMatrix[i*graphSize + j] == 1)      // there is a connection to that node 
				if (graphColors[j] == nodeColor) 
				{ 
					cout << "Color collision from node: " << i << " colored with: " << nodeColor << "  to node: " << j << " colored with " << graphColors[j] << endl; 
					numErrors++; 
					numErrorsOnRow++; 
				} 
		} 
		
		if (numErrorsOnRow != 0) 
			cout << "Errors for node " << i << " : " << numErrorsOnRow << endl; 
	} 
	
	cout << "Number of Color used: " << maxColor << " ; errors: " << numErrors << endl << "==================== " << endl ;    
} 



// Author: Pascal 
// Checks if a graph has conflicts or not from adjacency List
int checkCorrectColoring(int *adjacencyList, int *graphColors, int graphSize, int maxDegree, int interactive){
    int numErrors = 0;
    int maxColor = -1;	/* 最大颜色即总颜色数 */
	
	if (interactive != 2)
    	cout << endl << "==================" << endl << "Error checking for Graph" << endl;
	
    for (int i=0; i<graphSize; i++)                 // we check each row
    {
        int nodeColor = graphColors[i];
        int numErrorsOnRow = 0;

	if (nodeColor > maxColor)
		maxColor = nodeColor;
		
        for (int j=0; j<maxDegree;j++){
			
            if (adjacencyList[i*maxDegree + j] == -1)
                break;
            else{     // there is a connection to that node
                int node = adjacencyList[i*maxDegree + j];
                if (graphColors[node] == nodeColor)
                {
                    cout << "Color collision from node: " << i << " col with " << nodeColor << "    to: " << node << " col with " << graphColors[node] << endl;
                    numErrors++;
                    numErrorsOnRow++;
                }
            }
        }
		
        if (numErrorsOnRow != 0)
            cout << "Errors for node " << i << " : " << numErrorsOnRow << endl;
    }
	
	if (interactive != 2)
    	cout << "Number of colors used: " << maxColor << "  ; Color errors for graph : " << numErrors << endl << "==================== " << endl ;   
   
    return maxColor;
}




//----------------------- Other -----------------------//
// Any additional stuff needed ....
void getGraphStats(int *adjacencyList,  int graphSize, int maxDegree, int *startPartitionList, int *endPartitionList){
	int boundaryCount = 0;
	int boundaryIndex = 0;
	int partitionIndex = 0;
	int start, end, node;
	int hasOutside = 0;
	int countDeg;
	
	int *numInside, *numOutside, *degreeList;
	numInside = new int[graphSize];
	numOutside = new int[graphSize];
	degreeList = new int[graphSize];
	memset(numInside,0,graphSize*sizeof(int));
	memset(numOutside,0,graphSize*sizeof(int));
	memset(degreeList,0,graphSize*sizeof(int));
	
	for (int i=0; i<graphSize; i++){
		if (!(i < endPartitionList[partitionIndex]))
			partitionIndex++;
		
		start = startPartitionList[partitionIndex];
		end = endPartitionList[partitionIndex];
		
		countDeg = 0;
		for (int j=0; j<maxDegree; j++)
		{
			node = adjacencyList[i*maxDegree + j];
			if (node == -1)
				break;
			else{
				countDeg++;
				
				if ((node < start) || (node >= end))
					numOutside[i]++;
				else
					numInside[i]++;
			}
		}
		degreeList[i] = countDeg;
	}
	
	
	
	ofstream numInsideF ("numInside.csv");
	ofstream numOutsideF ("numOutside.csv");
	ofstream degreeListF ("degreeList.csv");
	
	for (int i=0; i<graphSize; i++){
		numInsideF << numInside[i] << " , ";
		numOutsideF << numOutside[i] << " , ";
		degreeListF << degreeList[i] << " , ";
	}
	
	numInsideF.close();
	numOutsideF.close();
	degreeListF.close();
	
	delete []numInside;
	delete []numOutside;
	delete []degreeList;
}

//----------------------- The meat -----------------------//

int main(int argc, char *argv[]){  
	
	if (! (((argc == 8) || (argc == 13)) || (argc == 16)) ){
	//if (((argc != 8) || (argc!=13))|| (argc !=16)){
		cout << "Arguments passed: " << argc << endl;
		cout << "8 or 13 or 16 Arguments needed:" << endl;
		cout << "gc <passes(0:auto/1:stated)> <input (0:real/1:artif)>  <metis (1:metis)>  <GPU rand(0-2)> <CPU algo (0:FF/1:SDO/2:part/3:tree/4:SDO_Out)> <GPU algo (0FF/1:SDO/2:Max/3:Min) inter(0:batch/1:interactive/2:csv-batch)" << endl;
		cout << "additional if batch: gridsize blocksize path+graphName weighted(y/n) less-than-deg-ok(y/n) " << endl;
		cout << "additional if batch and metis: metisInputFile numMetisPartitions metisOutput-to-use " << endl;
		
		return 1;
	}
	
	int maxDegree, numColorsSeq, numColorsParallel, boundaryCount, conflictCount, passes, graphSize, graphSizeRead;
	int _gridSize, _blockSize, numMetisPartitions, randomnessValue, avgDegree;
	float density;
	char ans;
	char weighted='i';
	long numEdges;
	string inputFilename;
	
	
	int *adjacentList, *adjacencyMatrix;

	conflictCount = boundaryCount = numColorsSeq = numColorsParallel = 0; 
	
	
	
	//--------------------- Parameter initialization ---------------------!
	bool useMetis = true;
	bool sdo = true;
	bool sdoConflictSolver = true;
	bool CPUSDO;
	int GPUSDO;
	int loopGPUS = 5;
	int GPUIterations = 10;		/* 并行着色迭代10次 */
	
	long randSeed = time(NULL);	
	randSeed = 1272167817;			// to set to a specific random seed for replicability
	
	
	passes = atoi(argv[1]);				// get number of passes
	randomnessValue = atoi(argv[4]);	// get randomness
	
	
	
	if (atoi(argv[3]) == 1)
		useMetis = true;
	else
		useMetis = false;
	
	
	
	bool artificial = false;
	if (atoi(argv[2]) == 0)
		artificial = false;
	else
		artificial = true;
	
	int seqTech;	
	seqTech = atoi(argv[5]);
	/*
	if (atoi(argv[5]) == 1)
		CPUSDO = true;
	else
		CPUSDO = false;
	*/
	
	GPUSDO = atoi(argv[6]);
	
	int interactive = atoi(argv[7]);

	
	//--------------------- Graph Creation ---------------------!
	
	// Grid and block size
	
	if (interactive != 2){
		cout << endl << "!--------------- Graph Coloring program -------------------!" << endl;
		cout << "Enter grid size (e.g 4): ";
	}
	if (interactive == 1)
		cin >> _gridSize;
	else{	/* interactive = 0 */
		_gridSize = atoi(argv[8]);

		if (interactive != 2)
			cout << _gridSize << endl;
	}

	if (interactive != 2)
		cout << "Enter block size (e.g 64): ";
	
	if (interactive == 1)
		cin >> _blockSize;
	else{
		_blockSize = atoi(argv[9]);
		
		if (interactive != 2)
			cout << _blockSize << endl;
	}
	if (interactive == 1){
		cout << endl << "Number of threads: " << _gridSize*_blockSize << endl;
		cout << endl;
	}



	
	int *startPartitionList = new int[_gridSize*_blockSize];	
	int *endPartitionList = new int[_gridSize*_blockSize];		
	memset(startPartitionList, -1, _gridSize*_blockSize*sizeof(int));
	memset(endPartitionList, -1, _gridSize*_blockSize*sizeof(int));
	
	int numRandoms = _gridSize*_blockSize*10;
	
	
	// Artificial or real 
	
	/* 使用真实graph */
	if (artificial == false){	// input file required - returns an adjacency list of the graph, graph size and max degree
		ifstream testFile;	
		
		if (interactive == 1){
			cout << "Enter graph input filename (e.g. graphs/1138_bus.mtx): ";
		//if (interactive == 1)
			cin >> inputFilename;
		}
		else
			inputFilename = argv[10];
		
		
		if ((interactive == 0)||(interactive == 2))
			weighted = argv[11][0];
		

		/*
		 // gets a compact adjacency list from the file input
		 readGraph(adjacencyMatrix, inputFilename.c_str(), _gridSize, _blockSize, graphSizeRead, graphSize, numEdges);
		 cout << graphSizeRead << " - " << graphSize << " - " << numEdges << endl; 
		 
		 
		 // gets the max degree
		 maxDegree = getMaxDegree(adjacencyMatrix, graphSize, avgDegree);
		 cout << "Max degree: " << maxDegree << "   average degree: " << avgDegree << endl;
		
		 
		 // Get adjacency list
		 adjacentList = new int[graphSize*maxDegree];
		 memset(adjacentList, -1, graphSize*maxDegree*sizeof(int)); 
		 
		 getAdjacentList(adjacencyMatrix, adjacentList, graphSize, maxDegree);
		 
		 
		 delete []adjacencyMatrix;
		 adjacencyMatrix = NULL;
		 */
		
		readGraphAdj(adjacentList, inputFilename.c_str(), _gridSize, _blockSize, graphSizeRead, graphSize, numEdges, maxDegree, weighted, interactive, useMetis);
		
		/* 输出邻接表到文件 */
		displayAdjacencyList(adjacentList, graphSize, maxDegree);
	}
	else
	{
		cout << "Enter graph size: ";
		cin >> graphSize;
		
		cout << "Enter density: ";
		cin >> density;
		
		numEdges = density*graphSize*(graphSize-1)/2;
		
		
		adjacencyMatrix = new int[graphSize*graphSize];  
		memset(adjacencyMatrix, 0, graphSize*graphSize*sizeof(int)); 
		
		
		srand ( randSeed );  // initialize random numbers  
		
		
		// generates a graph
		generateMatrix(adjacencyMatrix, graphSize, numEdges);
		
		
		// gets the max degree
		maxDegree = getMaxDegree(adjacencyMatrix, graphSize, avgDegree);
		cout << "Max degree: " << maxDegree << "   average degree: " << avgDegree << endl;
		
		
		// Get adjacency list
		adjacentList = new int[graphSize*maxDegree];
		memset(adjacentList, -1, graphSize*maxDegree*sizeof(int)); 
		
		getAdjacentList(adjacencyMatrix, adjacentList, graphSize, maxDegree);
		
		delete []adjacencyMatrix;
		adjacencyMatrix = NULL;
	}
	
	if (maxDegree > TEMP_COLOR_LENGTH)
	{
		if (interactive == 1){
			cout << endl << "Warning ... degree of graph (" << maxDegree << ") exceeds current TEMP_COLOR_LENGTH (" <<  TEMP_COLOR_LENGTH << ") " << endl;
			cout << "This might cause errors!!!" << endl;
			cout << "Are you sure you want to continue (y/n): ";
		//if (interactive == 1)
			cin >> ans;
		}
		else
			ans = argv[12][0];
		
			if (ans != 'y'){
				cout << "Exiting now - please change value of constant TEMP_COLOR_LENGTH in graphColoring.hi or use OPTION TWO in graphColoring.cu" << endl << endl;
				exit(0);
			}
		
	}
	else
		cout << "Allocation successful!" << endl;
	
	
	
	
	// Some further intializations
	int *graphColors = new int[graphSize];          
	int *boundaryList = new int[graphSize]; 
	
	//int *degreeList = new int[graphSize];
	degreeList = new int[graphSize];
	
	int *numOutside = new int[graphSize];
	
	memset(graphColors, 0, graphSize*sizeof(int)); 
	memset(boundaryList, 0, graphSize*sizeof(int)); 
	memset(degreeList, 0, graphSize*sizeof(int)); 
	memset(numOutside, 0, graphSize*sizeof(int));
	
	// Get degree List
	avgDegree = getDegreeList(adjacentList, degreeList, graphSize, maxDegree);
	
	
	
	
	int *randomList = new int[numRandoms];
	for (int i=0; i<numRandoms; i++)					// stores random numbers in the range of 0 to 2
		randomList[i] = rand()%(randomnessValue+1);
	
	
	//--------------------- Metis ---------------------!
	
	int *metisAdjacencyList = new int[graphSize*maxDegree];
	int *metisDegreeList = new int[graphSize];
	
	
	
	if (useMetis == true){
		memset(metisAdjacencyList, 0, graphSize*maxDegree*sizeof(int)); 
		memset(metisDegreeList, 0, graphSize*sizeof(int)); 
		
		string metisInputFile, metisOutputFile;
		int numMetisPartitions = _gridSize*_blockSize;
	
		if (interactive != 1){
			metisInputFile = argv[13];
			numMetisPartitions = atoi(argv[14]);
			metisOutputFile = argv[15];
		}
		
		numMetisPartitions = metis(adjacentList, metisAdjacencyList, graphSize, numEdges, maxDegree, startPartitionList, endPartitionList, numMetisPartitions, metisInputFile, metisOutputFile, interactive);	// Metis
		
		getDegreeList(metisAdjacencyList, metisDegreeList, graphSize, maxDegree);
		
		memcpy(adjacentList, metisAdjacencyList, graphSize*maxDegree*sizeof(int));
		memcpy(degreeList, metisDegreeList, graphSize*sizeof(int));
	}else{	/* 不使用metis */
		// allocating exact partition size for
		int partitionSizes = graphSize / (_gridSize*_blockSize);
		
		for (int i=0; i<_gridSize*_blockSize; i++){
			startPartitionList[i] = i*partitionSizes; 
			endPartitionList[i] = (i+1)*partitionSizes;
		}
	}
	
	
	
	
	//--------------------- Boundary List ---------------------!
	cudaEvent_t start_bc, stop_bc;
        float elapsedTimeBoundaryc = 0;

	//if (CPUSDO == true){
    if (seqTech != 0){
		cudaEventCreate(&start_bc);
        cudaEventCreate(&stop_bc);
        cudaEventRecord(start_bc, 0);

        getBoundaryNodes(adjacentList, boundaryList, graphSize, maxDegree, _gridSize,_blockSize, startPartitionList, endPartitionList);	

        cudaEventRecord(stop_bc, 0);
        cudaEventSynchronize(stop_bc);
        cudaEventElapsedTime(&elapsedTimeBoundaryc, start_bc, stop_bc);	
	}




	///////// Needs to be updated to use adjacencyList instead of adjacencyMatrix!!!!!!!!!!!!!!!!!!!
	cudaEvent_t start_b, stop_b;
	float elapsedTimeBoundary;
	cudaEventCreate(&start_b); 
	cudaEventCreate(&stop_b); 
	cudaEventRecord(start_b, 0); 
	
	//maxDegree = getBoundaryList(adjacencyMatrix, boundaryList, graphSize, boundaryCount, graphSize, _gridSize, _blockSize);	// return maxDegree + boundaryCount (as ref param)
	//boundaryCount = getBoundaryList(adjacentList, boundaryList, graphSize, maxDegree, _gridSize,_blockSize, startPartitionList, endPartitionList);				// get boundaryCount and get boundary list
	
	/* 计算边界结点 */
	if ((GPUSDO == 0) || (GPUSDO == 1))	/* FF || SDO 调用8个参赛的getBoundaryNodes */
		boundaryCount = getBoundaryNodes(adjacentList, boundaryList, graphSize, maxDegree, _gridSize,_blockSize, startPartitionList, endPartitionList);
	else	/* OMAX || OMIN 调用9个参赛的getBoundaryNodes */
		boundaryCount = getBoundaryNodes(adjacentList, boundaryList, numOutside, graphSize, maxDegree, _gridSize,_blockSize, startPartitionList, endPartitionList);				// get boundaryCount and get boundary list
	
	
	/*
	 cout << "Num outside: " << endl;
	 for (int i=0; i<graphSize; i++){
	 cout << i << " : " << numOutside[i] << endl; 
	 }
	 */
	
	
	
	cudaEventRecord(stop_b, 0); 
	cudaEventSynchronize(stop_b); 
	cudaEventElapsedTime(&elapsedTimeBoundary, start_b, stop_b); 	/* 计算边界结点的时间 */
	//cout << "Time to getBoundaryList :"<< elapsedTimeBoundary << " ms" << endl;
	
	
	//	for (int i=0; i<graphSize; i++)
	//		cout << i << ": " << boundaryList[i] << endl; 
	
	
	//--------------------- print graph statistics into a file ---------------------!
	getGraphStats(adjacentList, graphSize, maxDegree, startPartitionList, endPartitionList);
	
	
	
	//--------------------- Sequential Graph Coloring ---------------------!
	cudaEvent_t start, stop, stop_1, stop_4;         
	
	float elapsedTimeCPU, elapsedTimeGPU_1, elapsedTimeGPU_4; 
	float elapsedTimesGPU_1[MAXGPUITERATIONS], elapsedTimesGPU_4[MAXGPUITERATIONS]; 
	float elapsedTimesGPU[MAXGPUITERATIONS], elapsedTimesGPU_Bound[MAXGPUITERATIONS], elapsedTimesGPU_Col[MAXGPUITERATIONS], elapsedTimesGPU_Sol[MAXGPUITERATIONS]; 
	float elapsedTimeGPU, elapsedTimeGPU_Bound, elapsedTimeGPU_Col, elapsedTimeGPU_Sol; 
	int numColorsParallels[MAXGPUITERATIONS], conflictCounts[MAXGPUITERATIONS], interColorsParallels[MAXGPUITERATIONS];
	int avgConflicts, avgInterColorsParallel;
    float avgColors;	
	avgConflicts = avgColors = avgInterColorsParallel = 0;
	elapsedTimeGPU = 0;
	elapsedTimeGPU_Col = 0;
	elapsedTimeGPU_Bound = 0;
	elapsedTimeGPU_Sol = 0;

	
	
	cudaEventCreate(&start); 
	cudaEventCreate(&stop); 
	cudaEventRecord(start, 0);  
	
	
	// Original adjacency List
	/*
	if (CPUSDO == true)
		if (interactive != 2)
			sdoIm(adjacentList, graphColors, degreeList, graphSize, maxDegree);
			//seqGraphColoring(adjacentList, graphColors, degreeList, maxDegree, graphSize);
			//cout << "Iterations:" << sdoIm(adjacentList, graphColors, degreeList, graphSize, maxDegree,_gridSize*_blockSize, startPartitionList, endPartitionList) << endl;
		else
			//sdoIm(adjacentList, graphColors, degreeList, graphSize, maxDegree,_gridSize*_blockSize, startPartitionList, endPartitionList);
			seqGraphColoring(adjacentList, graphColors, degreeList,  maxDegree, graphSize);
	else
		colorGraph_FF(adjacentList, graphColors, graphSize, maxDegree);  
	*/

	if (seqTech == 1)
		sdoIm(adjacentList, graphColors, degreeList, graphSize, maxDegree);
	else
		if (seqTech == 2)
			sdoIm(adjacentList, graphColors, degreeList, graphSize, maxDegree,_gridSize*_blockSize, startPartitionList, endPartitionList);
		else
			if (seqTech == 3)
				seqGraphColoring(adjacentList, graphColors, degreeList, maxDegree, graphSize);
			else
				if (seqTech == 4)
					sdoImO(adjacentList, graphColors, degreeList, graphSize, maxDegree);
				else
					if (seqTech == 0)	/* 串行FF */
						colorGraph_FF(adjacentList, graphColors, graphSize, maxDegree);
    

	
	
	cudaEventRecord(stop, 0); 
	cudaEventSynchronize(stop); 
	cudaEventElapsedTime(&elapsedTimeCPU, start, stop); 	/* 串行着色的时间 */
	
	
//exit(0);	
	
	//--------------------- Checking for color conflict ---------------------!
	if (interactive != 2)
		cout << endl << endl << "Sequential Conflict check: ";
	
	/* 检查着色是否正确 */
	numColorsSeq = checkCorrectColoring(adjacentList, graphColors, graphSize, maxDegree,interactive);
	cout << "Sequential time: " << elapsedTimeCPU << " ms." << endl;

	if (interactive == 1){
		char ansYN;
		cout << "Do you want to stop now?: (y/n) ";
		cin >> ansYN;

		if (ansYN == 'y')
			exit(0);
	}
	
//	exit(0);
	
	
	//-------------------- Parallel Graph Coloring ---------------------!	
	
	int *conflict = new int[boundaryCount];                    // conflict array
	
	
	float GPU_Time = 0.0;
	float elapsedTimesCsFF = 0.0;
	/* 重复10次取平均值 */
	for (int repeatGPU=0; repeatGPU<GPUIterations; repeatGPU++)
	{
		memset(conflict, 0, boundaryCount*sizeof(int));            // conflict array initialized to 0  
		memset(graphColors, 0, graphSize*sizeof(int));             // reset colors to 0 
		
		//--------------- Steps 1, 2 & 3: Parallel Partitioning + Graph coloring + Conflict Detection
		cudaEventCreate(&start); 
		cudaEventCreate(&stop); 
		cudaEventCreate(&stop_1); 
		cudaEventCreate(&stop_4); 
		cudaEventRecord(start, 0); 
		
		
		int *conflictTmp = new int[boundaryCount];
		memset(conflictTmp, 0, boundaryCount*sizeof(int));  
		
		/* 调用kernel */
		GPU_Time += cudaGraphColoring(adjacentList, boundaryList, graphColors, degreeList, conflictTmp, boundaryCount, 
						  maxDegree, graphSize, passes, _gridSize*_blockSize, _gridSize, _blockSize,
						  startPartitionList, endPartitionList, randomList, numRandoms, GPUSDO, numOutside);
		
		cudaEventRecord(stop_1, 0); 
		cudaEventSynchronize(stop_1); 	/* stop_1 - start = GPU上的时间 */
		
		
		
		
		// count number of parallel colors
		/* 计算并行着色（未解决冲突）的颜色数 */
		int interColorsParallel = 0;
		for (int i=0; i<graphSize; i++)
			if ( interColorsParallel < graphColors[i] )
				interColorsParallel = graphColors[i];
		
		
		interColorsParallels[repeatGPU] = interColorsParallel;
		
		
		
		//-------- Conflict Count
		conflictCount = 0;
		for (int i=0; i<boundaryCount; i++)
		{
			int node = conflictTmp[i];
			
			if (node != -1)
			{
				//cout << "conflict " << conflictCount <<  " at " << node << endl;
				conflict[conflictCount] = node;
				conflictCount++;
			}
		}
		delete[] conflictTmp;
		
		
		conflictCounts[repeatGPU] = conflictCount;
		
		cudaEventRecord(stop_4, 0); 
		cudaEventSynchronize(stop_4); 	/* stop_4 - start = GPU上的时间 + 计算conflict结点的时间 */
		
		
		
		//--------------- Step 4: solve conflicts 
		
	
		cudaEvent_t csFF_start, csFF_stop;
		cudaEventCreate(&csFF_start); 
		cudaEventCreate(&csFF_stop);
		cudaEventRecord(csFF_start, 0); 
		/* CPU串行冲突解决 */
		if (GPUSDO == 0)	/* FF */
		    conflictSolveFF(adjacentList,  graphSize, conflict, conflictCount, graphColors, maxDegree);
		else	/* SDO */
			conflictSolveSDO(adjacentList, conflict, conflictCount, graphColors, degreeList, graphSize, maxDegree);
		cudaEventRecord(csFF_stop, 0); 
		cudaEventSynchronize(csFF_stop);
		
		cudaEventRecord(stop, 0); 
		cudaEventSynchronize(stop); 	/* stop - start = GPU上的时间 + 计算conflict结点的时间 + 串行冲突解决的时间 */
		float t_csFF;
		cudaEventElapsedTime(&t_csFF, csFF_start, csFF_stop);
		elapsedTimesCsFF += t_csFF;
		printf("csFF time: %f ms\n", t_csFF);
		
		cudaEventElapsedTime(&elapsedTimesGPU[repeatGPU], start, stop); 
		cudaEventElapsedTime(&elapsedTimesGPU_1[repeatGPU], start, stop_1); 
		cudaEventElapsedTime(&elapsedTimesGPU_4[repeatGPU], start, stop_4); 
		
		
		//--------------------- Checking for color conflict ---------------------!
		if (interactive != 2)
			cout << endl <<  "Parallel Conflict check:   Run: " << repeatGPU;	
		
		/* 检查着色是否正确并更新并行着色的颜色数 */
		numColorsParallels[repeatGPU] = checkCorrectColoring(adjacentList, graphColors, graphSize, maxDegree, interactive);
	}
	
	/* 各个时间分布 */
	int minColor = 100000;	
	for (int k=0; k<GPUIterations; k++){
		elapsedTimesGPU_Col[k] = elapsedTimesGPU_1[k];
		elapsedTimeGPU_Col += elapsedTimesGPU_Col[k];								// coloring time
		
		elapsedTimesGPU_Bound[k] = elapsedTimesGPU_4[k] - elapsedTimesGPU_1[k];		// getting number of conflicts
		elapsedTimeGPU_Bound += elapsedTimesGPU_Bound[k];
			
		elapsedTimesGPU_Sol[k] = elapsedTimesGPU[k] - elapsedTimesGPU_4[k];			// solving for conflicts
		elapsedTimeGPU_Sol += elapsedTimesGPU_Sol[k];
		
		elapsedTimeGPU += elapsedTimesGPU[k];										// total time
		
		
		avgInterColorsParallel += interColorsParallels[k];
		avgConflicts += conflictCounts[k];

		/* 最小颜色数 */
		if (minColor > numColorsParallels[k])
			minColor = numColorsParallels[k];

		avgColors += numColorsParallels[k];
	}
	
	
	/* 取平均值 */
	elapsedTimeGPU = elapsedTimeGPU/GPUIterations;
	elapsedTimeGPU_Col = elapsedTimeGPU_Col/GPUIterations; 
	elapsedTimeGPU_Bound = elapsedTimeGPU_Bound/GPUIterations;
	elapsedTimeGPU_Sol = elapsedTimeGPU_Sol/GPUIterations;
	avgConflicts = avgConflicts/GPUIterations;
	avgColors = avgColors/GPUIterations;
	avgInterColorsParallel = avgInterColorsParallel/GPUIterations;
	
	/*
	// Display average timings
	cout << endl << endl << "Coloring runs: " << endl;
	for (int k=0; k<GPUIterations; k++){
		cout << "others: " << elapsedTimesGPU_1[k] << " , " << elapsedTimesGPU_4[k]<< endl;   
		cout << "Run " << (k+1) << " : Coloring: " << elapsedTimesGPU_Col[k]  << " ms " << 
							       "    Boundary : " << elapsedTimesGPU_Bound[k]<< " ms " << 
								   "    Conflict : " << elapsedTimesGPU_Sol[k]  << " ms " << 
								   "    Total: "     << elapsedTimesGPU[k] + elapsedTimeBoundary  << " ms." << endl;	
		cout << "Conflicts: " <<  conflictCounts[k] << "   GPU Colors: " << interColorsParallels[k] << "   Colors: " << numColorsParallels[k] << endl << endl;
	}
	*/
	
//cout << graphSize << endl;	
	//--------------------- Information Output ---------------------!	
	if (interactive == 2){
		 float timeCPU;
        timeCPU = elapsedTimeCPU+elapsedTimeBoundaryc;

		if (GPUSDO == 0)
                	cout << "0 ";
		else
                        if (GPUSDO == 1)
			 	cout << "1 ";
                        else
                                if (GPUSDO == 2)
					cout << "2 ";
                                else
                                        if (GPUSDO == 3)
						cout << "3 ";

		cout << _gridSize << " " << _blockSize << " " << _gridSize*_blockSize << " " << graphSize/(_gridSize*_blockSize) << " " << passes << " " << avgConflicts << " " << timeCPU << " " << (elapsedTimeBoundary + elapsedTimeGPU) << " " << timeCPU/(elapsedTimeBoundary + elapsedTimeGPU)  << " " << minColor << " " << avgColors << " " << numColorsSeq  << " " << graphSize << ";" << endl; 
	}
	else{	/* 非交互，interactive=0 */
	cout << endl << endl << "!------- Graph Summary:" << endl;
	cout << "Vertices: " << graphSize << "   Edges: " << numEdges << "   Density: " << (2*numEdges)/((float)graphSize*(graphSize-1)) << endl;
	cout << "Max  Degree: " << maxDegree << "    Average degree: " << avgDegree << endl;
	if (artificial == false){
		cout << "Graph read in: " << inputFilename << endl;
		cout << "Vertices in graph: " << graphSizeRead << endl;
	}
	else
		cout << "Random seed used: " << randSeed << endl;
	cout << endl;
	
	
	cout << "Grid Size: " << _gridSize << "    Block Size: " << _blockSize << "     Total number of threads: " << _gridSize*_blockSize << endl;
	cout << "Graph average subsize: " << graphSize/(_gridSize*_blockSize) << endl;
	
	if (useMetis == true)
		cout << "Number of metis partitions: " << numMetisPartitions << endl;
	cout << endl;
	
	cout << "GPU Passes done: " << passes << endl;
	

	cout << endl << "Average timing for " <<  GPUIterations << " runs: " << endl;
	if (seqTech != 0){//CPUSDO == true){
		elapsedTimeCPU = elapsedTimeCPU+elapsedTimeBoundary;

		if (seqTech == 1)
			cout << "CPU time (SDO): " << elapsedTimeCPU << " ms ";
		else
			if (seqTech == 2)
				cout << "CPU time (SDO in parts): " << elapsedTimeCPU << " ms ";
			else
				if (seqTech == 4)
					cout << "CPU time (SDO Color allocation out): " << elapsedTimeCPU << " ms ";
				else
					cout << "CPU time (tree SDO): " << elapsedTimeCPU << " ms";

		if (GPUSDO == 0)
			cout << "   -  GPU Time (FF Solver): " << elapsedTimeGPU << " ms" << endl;
		else
			if (GPUSDO == 1)
				cout << "   -  GPU Time (SDO Solver): " << elapsedTimeGPU << " ms" << endl;
			else
				if (GPUSDO == 2)
					cout << "   -  GPU Time (Max SDO Solver): " << elapsedTimeGPU << " ms" << endl;
				else
					if (GPUSDO == 3)
						cout << "   -  GPU Time (Min SDO Solver): " << elapsedTimeGPU << " ms" << endl;
	}
	else	/* FF */
	{
		cout << "CPU time (FF): " << elapsedTimeCPU << " ms ";
		if (GPUSDO == 0)	/* GPU FF */
			cout << "   -  GPU Time (FF Solver): " << elapsedTimeGPU << " ms" << endl;
		else
			if (GPUSDO == 1)
				cout << "   -  GPU Time (SDO Solver): " << elapsedTimeGPU << " ms" << endl;
			else
				if (GPUSDO == 2)
					cout << "   -  GPU Time (Max SDO Solver): " << elapsedTimeGPU << " ms" << endl;
				else
					if (GPUSDO == 3)
						cout << "   -  GPU Time (Min SDO Solver): " << elapsedTimeGPU << " ms" << endl;
	}
	
	cout << endl << "Getting boundary list: " 	<< elapsedTimeBoundary << " ms" << endl; 
	cout << "ALGO step 1, 2 & 3 (GPU)      : " 	<< elapsedTimeGPU_Col << " ms" << endl;  
	cout << "Count Conflicts               : " 	<< elapsedTimeGPU_Bound << " ms" << endl; 
	cout << "ALGO step 4 (solve conflicts) : " 	<< elapsedTimeGPU_Sol << " ms" << endl; 
	cout << "Total time                    : "	<< (elapsedTimeBoundary + elapsedTimeGPU) << " ms" << endl;
	cout << endl;
	
	
	cout << "Boundary Count: " << boundaryCount << endl;
	cout << "Avg Conflict count: " << avgConflicts << endl;
	cout << endl;
	
	cout << "Sequential Color: " << numColorsSeq << endl;
	cout << "Avg Colors before solving conflict: " << avgInterColorsParallel << endl;
	cout << "Min color allocated: " << minColor << "   -   Average color: " << avgColors << endl;
	cout << "Avg GPU speed up (including boundary): " << elapsedTimeCPU/(elapsedTimeBoundary + elapsedTimeGPU) << " x" << endl;
	
	cout << "!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!" << endl << endl;
	
	}

	GPU_Time /= GPUIterations;
	elapsedTimesCsFF /= GPUIterations;
	float parallel_time_total;
	parallel_time_total = GPU_Time + elapsedTimeBoundary + elapsedTimesCsFF;
	printf("Paralle: GPU_Time=%f, elapsedTimeBoundary=%f, elapsedTimesCsFF=%f, parallel_time_total=%f\n", GPU_Time, elapsedTimeBoundary, elapsedTimesCsFF, parallel_time_total);
	
	//--------------------- Cleanup ---------------------!		
	delete []graphColors; 
	delete []conflict; 
	delete []boundaryList;
	delete []adjacentList;
	delete []degreeList;
	delete []randomList;
	
	return 0;  
}  

