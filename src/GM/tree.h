
class node{
private:
	int key, saturation, degree, color;
	node *left, *right;
	
public:
	node();
	node(int index, int sat, int deg);
	node(int index, int sat, int deg, int col, node *L, node *R);
	
	
	int getKey();
	int getSaturation();
	int getDegree();
	int getColor();
	node* getLeft();
	node* getRight();
	
	void setKey(int index);
	void setSaturation(int sat);
	void setDegree(int deg);
	void setColor(int c);
	void setLeft(node *L);
	void setRight(node *R);
	
	void setKSD(int index, int saturation, int degree);
	
	void displayNode();
	
	~node();
};


// Tree sorted by saturation and then  by degree
class tree{
private:
	node *top;
	
public:
	tree();
	
	void insert(node *x);
	node* remove(int index, int saturation, int degree);
	
	node* findNode(int index, int saturation, int degree);
	void findBiggest(int &index, int &saturation, int &degree);
	
	void displayTreeRML(node *current);
	void displayTreeLMR(node *current);
	void displayTreeMLR(node *current);
	
	node* getTop();
	
	~tree();
};

