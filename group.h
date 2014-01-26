#include <assert.h>
#include <fstream>
#include <iostream>
#include <string>
#include <list>

using namespace std;

const float HUGEREAL = 100000;
const double TOL = 0.00001;
class Edge;
class Triangle;
class EdgeList;
class TriList;
class Ray;
class Color;
template< class NodeClass > class genList;

inline float myabs(float a)
{
	if(a < 0)
		return -a;
	else return a;
}

class Vector
{
	friend ostream& operator<<(ostream&, const Vector&);	// printing in Vector form
public:
	Vector(){x=0;y=0;z=0;}
	Vector(float a, float b, float c) 
	{
		x = a;
		y = b;
		z = c;
	}
	float x, y, z;

	/* operator overloading functions */
	Vector operator*(const Vector&);	// Vector cross product
	Vector operator+(const Vector&);	// Vector addition
	Vector operator-(const Vector&);	// Vector subtraction
	Vector operator~();					// unit Vector
	Vector scale(float t);				// scalar product
	float operator&(const Vector&);		// Vector dot product
	float mod();						// magnitude of Vector0
};

class Color
{
public:
	Color(): red(0.0), green(0.0), blue(0.0) {}
	float red;
	float green;
	float blue;
};

class Vertex
{
	friend class VertTree;
public:
	Vertex(const Vector& vect):leftptr(0), rightptr(0), next(0), v(vect), trinum(1), edgeptr(0) {}
	Vector getVector() const { return v;}
	void setVector( const Vector& setv) { v = setv; } 
	void setEdge(Edge* newptr) {this->edgeptr = newptr;}
	Edge* getEdge() { return edgeptr; };
	void setAvgNorm( Vector norm) { avg_norm = norm; }
	Vector getAvgNorm() const { return avg_norm; }
private:
	Vertex* leftptr;	// pointer to left subtree
	Vertex* rightptr;	// pointer to right subtree
	Vertex* next;		// forward pointer to copy
	Vector v;			// vertex coordinates
	int trinum;			// # of incident triangles or # of copies
	Edge* edgeptr;		// pointer to half-edge emanating from the vertex
	Vector avg_norm;	// average normal (computed at run-time before firing rays)
};

class VertTree
{
public:
	VertTree() {rootptr = 0;};
	void insertNode(const Vector&, const int i, int&);
	void inOrderTraversal() const;
private:
	Vertex* rootptr;
	void insertNodeHelper(Vertex**, const Vector&, const int i, int&);
	void inOrderHelper(Vertex*) const;
};

class Edge
{
	friend class EdgeList;
public:
	Edge():nextptr(0),verta(0), vertb(0), next(0), pair(0), face(0) {}
	void insertEdge(Vertex*, const int, int&);
	Vertex* getVerta() const { return verta; }
	Vertex* getVertb() const { return vertb; }
	Edge* getNext() const { return next; }
	Edge* getPair() const { return pair; }
	void setFace(Triangle* facet) {face = facet;}
	Triangle* getFace() const { return face; }
	void setPair(Edge* edge) { pair = edge; }
	Edge* nextptr;
private:
	Vertex* verta;
	Vertex* vertb;
	Edge* next;
	Edge* pair;
	Triangle* face;
};

class Triangle
{
	friend class TriList;
public: 
	Triangle(): bits(0), nextptr(0), group_num(0) {}
	Vector Blend_norm(Vector);
	void setnormal(const Vector vect) { facet_norm = vect; }
	Vector getnormal() { return facet_norm; }
	void setVerts(Vertex** vert);
	Vertex** getVerts() { return vert; } 
	//Vertex* getVert() const { return vert[0]; } 
	void setEdges(Edge*, Edge*, Edge*);
	void setOneEdge(Edge* edge1) { edge = edge1; }
	Edge** getEdges() { return edges; } 
	int getGroupNum() const { return group_num; }
	Triangle* nextptr;	
	void setObjNum(const int num) { objnum = num; }
	const int getObjNum() const { return objnum; } 
private:
	int bits;			// # of times each triangle is visited in group# traversal
	int group_num;		// triangle group number
	int objnum;
	Vector facet_norm;  // facet normal
	Vertex* vert[3];	// vertices
	Edge* edges[3];		// 3 pointers to all edges 
	Edge* edge;			// pointer to first edge
};

class Object			// object description
{
public:
	Object(): nextptr(0),objname(""), firstTri(0), lastTri(0) {}
	void setbound(Vector max, Vector min)
	{
		vmax = max;
		vmin = min;
	}
	float bound_box(Ray) const;
	float Shadow_box(Ray ray) const;
	Vector getMinVert() { return vmin; }
	Vector getMaxVert() { return vmax; }
	Object* nextptr;
	void setFirstTri(Triangle* firstptr) { firstTri = firstptr; }
	void setLastTri(Triangle* lastptr) { lastTri = lastptr; }
	Triangle* getFirstTri() { return firstTri; }
	Triangle* getLastTri() { return lastTri; }
	void setname(string name) { objname = name; }
	string getname() { return objname; }
	void setObjNum(const int num) { objnum = num; }
	const int getObjNum() { return objnum; }
	void setIamb(Color Iamb) { iamb = Iamb; }
	Color getIamb() { return iamb; }
	void setKdiff(Color Kdiff) { Kd = Kdiff; }
	Color getKdiff() { return Kd; }
	void setKspec(float Kspec) { Ks = Kspec; }
	float getKspec() { return Ks; }
	void setKref(Color Kref) { Kr = Kref; }
	Color getKref() { return Kr; }
	void setKtrans(Color Ktrans) { Kt = Ktrans; }
	Color getKtrans() { return Kt; }
	void setIndex(float N) { index = N; }
	float getIndex() { return index; }
	void setN(float N) { n = N; }
	float getN() { return n; }
	Vector vmin, vmax;	// for bounding box
private:
	string objname;		// object name
	Triangle* firstTri;	// pointer to first triangle
	Triangle* lastTri;	// pointer to llast triangle
	//Vector vmin, vmax;	// for bounding box
	int objnum;
	Color iamb, Kd, Kr, Kt;
	float n;
	float Ks;
	float index;
};

template< class NodeClass >
class genList
{
public:
	genList();	
	//~genList() {};					// Destructor
	~genList();
	void insertnode(NodeClass* newptr);
	bool isEmpty() const;
	void print() const;
//private:
	NodeClass* firstptr;
	NodeClass* lastptr;
	
};

template< class NodeClass >	// default constructor
genList< NodeClass >::genList():firstptr(0), lastptr(0) {}

template< class NodeClass >
genList<NodeClass>::~genList()
{
	if(!isEmpty())
	{
		NodeClass* currentptr = firstptr;
		NodeClass* tempptr;
		
		while(currentptr != 0)
		{
			tempptr = currentptr;
			currentptr = currentptr->nextptr;
			delete tempptr;
		}
	}
	
}

template< class NodeClass >
bool genList<NodeClass>::isEmpty() const
{ return firstptr == 0; }

template< class NodeClass >
void genList<NodeClass>::insertnode(NodeClass* newptr) 
{
	if(isEmpty())
		firstptr = lastptr = newptr;
	else
	{
		lastptr->nextptr = newptr;
		lastptr = newptr;
	}
}

template<class NodeClass>
void genList<NodeClass>::print() const
{
	if(isEmpty())
	{
		cout<<"The list is empty.\n";
		return;
	}
	
	NodeClass* currentptr = firstptr;
	
	while(currentptr != 0)
	{
		cout<<currentptr<<endl;
		currentptr = currentptr->nextptr;
	}
}

class EdgeList
{
public:
	EdgeList();	
	~EdgeList() {};	// Destructor
	
	void insertnode(Edge* newptr);
	bool isEmpty() const;
	void print() const;
	void findPair();
//private:
	Edge* firstptr;
	Edge* lastptr;
};

class TriList
{
public:
	TriList();	
	//~TriList() {};	// Destructor
	~TriList();
	void insertnode(Triangle* newptr);
	void AvgNorm() const;
	void computeGroup(Triangle* );
	bool isEmpty() const;
	void print() const;
//private:
	Triangle* firstptr;
	Triangle* lastptr;
	
};

class Ray
{
public:
	Ray(): depth(0), object(0) {} 
	float intersect_plane(Vector, Vector);
	int Compute_intersect();
	Ray Reflect_ray();
	Ray Transmit_ray();
	Color Shadow_ray();
	void setBase(const Vector eye) { base = eye; }
	void setDir(const Vector aim) { dir = aim; }
	void setObject(Object* obj) { object = obj; }
	int getDepth() { return depth; }
	void setDepth(int d) { depth = d; }
	Vector getBase() const { return base; }
	Vector getDir() const { return dir; }
	Object* getObject() const { return object; }
	
private:
	
	Vector base;		// Origin of ray
	Vector dir;			// Ray direction
	int depth;			// Level of recursion or depth of tree
	Vector hitpoint;
	Triangle* triangle;// Point of intersection
	Object* object;		// Object hit by ray
	Vector blend_norm;	// Blended normal of triangle hit by ray

};

struct ObjCopy
{
	float tvalue;
	Object* obj;
};

struct Light
{
	Vector v;
	Color col;
};