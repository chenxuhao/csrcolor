#include "tree.h"
#include <iostream>

using namespace std;


/*
int main(){
	tree graph;
	node *temp;
	
	
	node *nodes = new node[7];
	int degList[7] = {4,2,4,5,10,4,7};
	
	for (int i=0; i<7; i++){
		nodes[i].setKSD(i, 0, degList[i]);
		graph.insert(&nodes[i]);
	}
	
	
	cout << "RML" << endl;
	graph.displayTreeRML(graph.getTop()); 
	
	
	
	temp = graph.remove(1,0,2);
	cout << endl << "Deleted: "; temp->displayNode();
	
	cout << "RML" << endl;
	graph.displayTreeRML(graph.getTop());
	
	
	cout << endl << "Deleted: "; temp = graph.remove(6,0,7);
	temp->displayNode();
	
	cout << "RML" << endl;
	graph.displayTreeRML(graph.getTop());
	
	
	return 0;
}
*/



node::node(){
	key = saturation = degree = color = -1;
	
	left = NULL;
	right = NULL;
}

node::node(int index, int sat, int deg){
	key = index;
	saturation = sat;
	degree = deg;
	
	color = -1;
	
	left = NULL;
	right = NULL;
}

node::node(int index, int sat, int deg, int col, node *L, node *R){
	key = index;
	saturation = sat;
	degree = deg;
	color = col;
	left = L;
	right = R;
}


int node::getKey(){
	return key;
}

int node::getSaturation(){
	return saturation;
}

int node::getDegree(){
	return degree;
}

int node::getColor(){
	return color;
}

node* node::getLeft(){
	return left;
}

node* node::getRight(){
	return right;
}


void node::setKey(int index){
	key = index;
}

void node::setSaturation(int sat){
	saturation = sat;
}

void node::setDegree(int deg){
	degree = deg;
}

void node::setKSD(int index, int sat, int deg){
	key = index;
	saturation = sat;
	degree = deg;
	color = -1;
	left = NULL;
	right = NULL;
}

void node::setColor(int c){
	color = c;
}

void node::setLeft(node *L){
	left = L;
}

void node::setRight(node *R){
	right = R;
}


void node::displayNode(){
	cout <<  key  << " :  Sat: " << saturation << " , Deg: " << degree << endl; 
}


node::~node(){
	left = NULL;
	right = NULL;
}








tree::tree(){
	top = NULL;
}

void tree::insert(node *x){
	node *current, *previous;
	bool left;
	
	if (top == NULL)
		top = x;
	else 
	{
		current = top;

		// Check to see where to insert
		while (current != NULL){
			previous = current;
			
			if (current->getSaturation() < x->getSaturation()){				
				current = current->getRight();
				left = false;
			}
			else 
				if (current->getSaturation() > x->getSaturation()){
					current = current->getLeft();
					left = true;
				}
				else 
					if (current->getDegree() < x->getDegree()){
						current = current->getRight();
						left = false;
					}
					else 
						if (current->getDegree() >= x->getDegree()){
							current = current->getLeft();
							left = true;
						}
		}
		
		// Insert item
		if (left == true)
			previous->setLeft(x);
		else 
			previous->setRight(x);
	}
}


node* tree::findNode(int index, int saturation, int degree){
	node *current, *previous;
	
	current = top;
	
	while ((current != NULL) && (current->getKey() != index)){
		previous = current;
		
		if (current->getSaturation() < saturation)	
			current = current->getRight();
		else 
			if (current->getSaturation() > saturation)
				current = current->getLeft();
			else 
				if (current->getDegree() < degree)
					current = current->getRight();
				else 
					if (current->getDegree() >= degree)
						current = current->getLeft();
	}
	
	return current;
}

node* tree::remove(int index, int saturation, int degree){
	node *current, *previous, *parent, *nodeToDel;
	bool left, parentLeft;
	parent = previous = current = top;
	//node blank;
	
	if (top == NULL){
		cout << "Tree is empty!!!" << endl;
	}
	else{
		// step1: find the node
		while ((current != NULL) && (current->getKey() != index)){
			previous = current;
			
			if (current->getSaturation() < saturation){				
				current = current->getRight();
				left = false;
			}
			else 
				if (current->getSaturation() > saturation){
					current = current->getLeft();
					left = true;
				}
				else 
					if (current->getDegree() < degree){
						current = current->getRight();
						left = false;
					}
					else 
						if (current->getDegree() >= degree){
							current = current->getLeft();
							left = true;
						}
		}
		
		
		// Not found!!!
		if (current == NULL){
			cout << "Node not found!!!" << endl;
			return NULL;
		}
		
		
		
		// Replace
		parent = previous;
		parentLeft = left;
		nodeToDel = current;
		
	//	blank.setKey(nodeToDel->getKey());
	//	blank.setSaturation(nodeToDel->getSaturation());
	//	blank.setDegree(nodeToDel->getDegree());
		
		
		// Option 1: A leaf; replace by nothing!!!
		if ((current->getLeft() == NULL) && (current->getRight() == NULL)){
			if (parentLeft == true)
				parent->setLeft(NULL);
			else 
				parent->setRight(NULL);
			
			if (top == nodeToDel)
				top = NULL;
			
			return nodeToDel;
		}
		
		//Option 2: Node had only 1 child
		if ((current->getLeft() == NULL) || (current->getRight() == NULL)){
			if (current->getLeft() == NULL) 
				current = current->getRight();
			else 
				current = current->getLeft();
			
			if (top == nodeToDel)
				top = current;
			else
				if (parentLeft == true)
					parent->setLeft(current);
				else 
					parent->setRight(current);
			
			return nodeToDel;
		}
		
		
		
		//Option 3: Node had 2 Children - the painful one: replace by node slightly biggest (normally the rightmost of the left node)
		previous = current;
		current = current->getLeft();
		
		if (current->getRight() == NULL){
			if (top == nodeToDel)
				top = current;
			else
				if (parentLeft == true)
					parent->setLeft(current);
				else 
					parent->setRight(current);
			
			current->setRight(nodeToDel->getRight());
		}
		else{
			while (current->getRight() != NULL){
				previous = current;
				current = current->getRight();
			}
			
			
			if (current->getLeft() == NULL)	// replaced node is a leaf
				previous->setRight(NULL);
			else 
				previous->setRight(current->getLeft());	// node has left children
			
			
			current->setLeft(nodeToDel->getLeft());
			current->setRight(nodeToDel->getRight());
			
			if (top == nodeToDel)
				top = current;
			else
				if (parentLeft == true)
					parent->setLeft(current);
				else 
					parent->setRight(current);
		}
		
		return nodeToDel;
	}
}


void tree::findBiggest(int &index, int &saturation, int &degree){
	node *current, temp;
	
	current = top;
	
	while (current->getRight() != NULL){
		current = current->getRight();
	}
	
	index = current->getKey();
	saturation = current->getSaturation();
	degree = current->getDegree();
}



void tree::displayTreeRML(node *current){
	
	if (current != NULL) 
	{
		displayTreeRML(current->getRight());
		current->displayNode();
		displayTreeRML(current->getLeft());
	}
}

void tree::displayTreeLMR(node *current){
	
	if (current != NULL) 
	{
		displayTreeLMR(current->getLeft());
		current->displayNode();
		displayTreeLMR(current->getRight());
	}
}

void tree::displayTreeMLR(node *current){
	
	if (current != NULL) 
	{
		current->displayNode();
		displayTreeMLR(current->getLeft());
		displayTreeMLR(current->getRight());
	}
}


node* tree::getTop(){
	return top;
}

tree::~tree(){
	top = NULL;
}

