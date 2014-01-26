#include "group.h"
#include <fstream>
#include <iostream>
#include <string>
#include <list>
using namespace std;

const int SIZE = 255;
extern EdgeList edgelist;
extern TriList trianglelist;
genList<Object> objectlist;
int triflag;
int objnum = 0;

Vector norm;
VertTree tree;
Triangle* tri;
Object* obj;
Vector vmin(HUGEREAL, HUGEREAL, HUGEREAL), vmax(-HUGEREAL, -HUGEREAL, -HUGEREAL);

void Readprop(const char*);
void readfile(const char*);
void Ray_trace();

void main()
{
	char* filename;
	//filename = "tetrahedron.stl";
	filename = "t1.stl";
	//filename = "view2.txt";
	//filename = "sphere3.txt";
	//filename = "knob.txt";

	readfile(filename);

	//trianglelist.print();
	//cout<<endl;
	//edgelist.print();
	
	edgelist.findPair();

	//char* fileprop = "viewprop2.txt";
	//char* fileprop = "properties.txt";
	char* fileprop = "teaprop.txt";
	Readprop(fileprop);

	int count=0;
	
	trianglelist.AvgNorm();

	Object* objcount;
	objcount = objectlist.firstptr;
	while(objcount != 0)
	{
		cout<<objcount<<endl;
		Triangle* tri = objcount->getFirstTri();
		int num = objcount->getObjNum();
		//cout<<num<<endl;
		//trianglelist.computeGroup(tri);	// Group# method
		cout<<objcount->getMinVert()<<" "<<objcount->getMaxVert()<<endl;
		objcount = objcount->nextptr;
		
	}
	
	/*Triangle* t1;
	t1 = trianglelist.firstptr;
	cout<<endl;
	int tricount = 0;
	while(t1 != 0)
	{
		tricount++;
		int num = t1->getGroupNum();
		//cout<<tricount<<" "<<num<<" "<<t1->getObjNum()<<endl;
		t1 = t1->nextptr;
	}*/

	//tree.inOrderTraversal();	// called for averaging normals using group #s

	Ray_trace();	// This one call does everything
}

void readfile(const char* filename)
{
	//Vertex* vert[3];
	Triangle *first, *last;
	Vector vect, newvect;
	string keyword;
	string objname;
	float nx, ny, nz;
	float vx, vy, vz;
	int i;

	char line[SIZE];
	
	ifstream input_file(filename, ios::in);
	if(!input_file.is_open())
	{
		cout <<"Error: Cannot input filename: "<< filename << endl;
		exit(0);
	}

	while(input_file.peek() != EOF)
	{
		int j = 0;
		char c;
		keyword = "";

		while(isalnum(c = input_file.get()))
		{
			line[j] = c;
			j++;
		}
		line[j] = NULL;
		keyword = line;
		//cout<<keyword<<endl;

		if((keyword == "SOLID") || (keyword == "solid"))
		{
			input_file.getline(line, SIZE, '\n');
			objname = line;
			obj = new Object;
			obj->setname(objname);
			objnum++;
			obj->setObjNum(objnum);
			triflag = 1;
			cout<<objname<<endl;
		}

		if((keyword == "FACET") || (keyword == "facet"))
		{
			tri = new Triangle;
			tri->setObjNum(objnum);
			input_file.getline(line, SIZE, ' ');
			input_file >> nx >> ny >> nz;
			vect.x = nx;
			vect.y = ny;
			vect.z = nz;
			//tri->setnormal(vect);	// Assigning facet normals to triangles from STL file
			//cout<<tri<<" "<<tri->getObjNum()<<endl;
			input_file.getline(line, SIZE, '\n');
		}

		if((keyword == "OUTER") || (keyword == "outer"))
		{
			input_file.getline(line, SIZE, '\n');
			i = 0;
		}

		if((keyword == "VERTEX") || (keyword == "vertex"))
		{
			Vector vtol(TOL, TOL, TOL);

			input_file >> vx >> vy >> vz;
			if(vx < vmin.x)
				vmin.x = vx;
			if(vy < vmin.y)
				vmin.y = vy;
			if(vz < vmin.z)
				vmin.z = vz;
			if(vx > vmax.x)
				vmax.x = vx;
			if(vy > vmax.y)
				vmax.y = vy;
			if(vz > vmax.z)
				vmax.z = vz;
			vect.x = vx;
			vect.y = vy;
			vect.z = vz;
			tree.insertNode(vect, i, triflag);
			
			i++;
			input_file.getline(line, SIZE, '\n');
		}
				

		if((keyword == "ENDLOOP") || (keyword == "endloop"))
		{
			//cout<<keyword<<" reached"<<endl;
			//input_file.getline(line, SIZE, '\n');
		}
			
		if((keyword == "ENDFACET") || (keyword == "endfacet"))
		{
			trianglelist.insertnode(tri);
			//cout<<keyword<<" reached"<<endl;
			//input_file.getline(line, SIZE, '\n');
		}

		if((keyword == "ENDSOLID") || (keyword == "endsolid"))
		{
			input_file.getline(line, SIZE, '\n');
			objname = line;
			
			obj->setLastTri(tri);
			obj->setbound(vmax, vmin);
			objectlist.insertnode(obj);

			vmax.x = -HUGEREAL;
			vmax.y = -HUGEREAL;
			vmax.z = -HUGEREAL;
			vmin.x = HUGEREAL;
			vmin.y = HUGEREAL;
			vmin.z = HUGEREAL;

			first = obj->getFirstTri();
			last = obj->getLastTri();
			cout<<first<<" "<<last<<" here"<<endl;
			
			cout<<objname<<" "<<obj<<endl;
		}
		
	}
}