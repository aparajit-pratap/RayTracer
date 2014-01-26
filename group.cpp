#include <assert.h>
#include "group.h"
#include <fstream>
#include <iostream>
#include <string>
#include <list>

using namespace std;

Vertex* vert[3];

EdgeList edgelist;
TriList trianglelist;
extern genList<Object> objectlist;
extern Object* obj;
extern Triangle* tri;
extern int triflag;

Edge* edge1;
Edge* edge2;
Edge* edge3;

void VertTree::insertNode(const Vector& vect, const int i, int& triflag)
{
	insertNodeHelper(&rootptr, vect, i, triflag);
}

void VertTree::inOrderTraversal() const	// called for grouping
{
	inOrderHelper(rootptr);
}

void VertTree::insertNodeHelper(Vertex **ptr, const Vector& vect, const int i, int& triflag)
{
	Edge newEdge;
	if((*ptr == 0))
	{
		*ptr = new Vertex(vect);
		//cout<<(*ptr)->v.x<<" "<<(*ptr)->v.y<<" "<<(*ptr)->v.z<<endl;
		assert(ptr != 0);
		newEdge.insertEdge(*ptr, i, triflag);
	}

	else if( ((*ptr)->v.x - vect.x) > TOL )
		insertNodeHelper( &((*ptr)->leftptr), vect, i, triflag);
				
	else if( (vect.x - (*ptr)->v.x) > TOL )
		insertNodeHelper( &((*ptr)->rightptr), vect, i, triflag);
			
	else if( ((*ptr)->v.y - vect.y) > TOL )
		insertNodeHelper( &((*ptr)->leftptr), vect, i, triflag);
			
	else if( (vect.y - (*ptr)->v.y) > TOL )
		insertNodeHelper( &((*ptr)->rightptr), vect, i, triflag);
			
	else if( ((*ptr)->v.z - vect.z) > TOL )
		insertNodeHelper( &((*ptr)->leftptr), vect, i, triflag);
			
	else if( (vect.z - (*ptr)->v.z) > TOL )
		insertNodeHelper( &((*ptr)->rightptr), vect, i, triflag);
				
	else
	{
		(*ptr)->trinum++;
		insertNodeHelper(&((*ptr)->next), vect, i, triflag);
	}
}

// Computes average normals using group #s
void VertTree::inOrderHelper(Vertex *ptr) const	// called by inOrderTraversal()
{
	float xx;
	if(ptr != 0)
	{
		Edge *edge, *pair;
		int group, firstgroup;
		int count = 1;
		Triangle *face;
		Vector sum;

		edge = ptr->edgeptr;
		face = edge->getFace();
		firstgroup = face->getGroupNum();
		sum = face->getnormal();
		pair = edge->getPair();
		edge = pair->getNext();
		face = edge->getFace();
	
		while(edge != ptr->edgeptr)
		{
			if(group == firstgroup)
			{
				sum = sum + face->getnormal();
				count++;
			}
			pair = edge->getPair();
			edge = pair->getNext();
			face = edge->getFace();
			group = face->getGroupNum();
		}
		ptr->avg_norm = sum.scale(1/(float) count);
		xx = ptr->avg_norm.mod();
		//cout<<ptr->v.x<<" "<<ptr->v.y<<" "<<ptr->v.z<<" "<<ptr->trinum<<" "<<endl;
		cout<<ptr->avg_norm.x<<" "<<ptr->avg_norm.y<<" "<<ptr->avg_norm.z<<endl;

		if(ptr->next != 0)
			inOrderHelper(ptr->next);
		inOrderHelper(ptr->leftptr);
		inOrderHelper(ptr->rightptr);
	}
}


void Edge::insertEdge(Vertex* vtx, const int i, int& triflag)
{
	Edge* tempEdge;
	Edge* newEdge = new Edge[3];

	switch(i)
	{
	case 0:{
		edge1 = &newEdge[0];
		vtx->setEdge(edge1);
		break;
		   }
	case 1:{
		edge2 = &newEdge[1];
		vtx->setEdge(edge2);
		break;
		   }
	case 2:{
		edge3 = &newEdge[2];
		vtx->setEdge(edge3);
		break;
		   }
	}
	
	Vector v1, v2, v3, norm;
	vert[i] = vtx;
		
	tempEdge = vtx->getEdge();
	
	if(i==2)
	{
		edge1->verta = vert[0];
		edge1->vertb = vert[1];
		edge2->verta = vert[1];
		edge2->vertb = vert[2];
		edge3->verta = vert[2];
		edge3->vertb = vert[0];

		edge1->next = edge2;
		edge2->next = edge3;
		edge3->next = edge1;

		edge1->setFace(tri);
		edge2->setFace(tri);
		edge3->setFace(tri);

		edgelist.insertnode(edge1);
		edgelist.insertnode(edge2);
		edgelist.insertnode(edge3);

		v1 = vert[0]->getVector();
		v2 = vert[1]->getVector();
		v3 = vert[2]->getVector();
		norm = ~((v2 - v1)*(v3 - v1));

		tri->setnormal(norm);
		tri->setOneEdge(edge1);
		tri->setEdges(edge1, edge2, edge3);
		tri->setVerts(vert);
		
		if(triflag == 1)
			obj->setFirstTri(tri);
		
		triflag = 0;
		//trianglelist.insertnode(tri);
		
		/*
		cout<<&newEdge[0]<<endl;
		cout<<&newEdge[1]<<endl;
		cout<<&newEdge[2]<<endl;


		cout<<newEdge[0].verta<<" "<<newEdge[0].vertb<<endl;
		cout<<newEdge[1].verta<<" "<<newEdge[1].vertb<<endl;
		cout<<newEdge[2].verta<<" "<<newEdge[2].vertb<<endl;*/
						
	}
}

void Triangle::setEdges(Edge* edge1, Edge* edge2, Edge* edge3)
{
	this->edges[0] = edge1;
	this->edges[1] = edge2;
	this->edges[2] = edge3;
}

void Triangle::setVerts(Vertex** vert)
{
	this->vert[0] = vert[0];
	this->vert[1] = vert[1];
	this->vert[2] = vert[2];
}


/*************************************************************************/

EdgeList::EdgeList():firstptr(0), lastptr(0) {}

/*EdgeList::~EdgeList()
{
	if(!isEmpty())
	{
		Edge* currentptr = firstptr;
		Edge* tempptr;
		
		while(currentptr != 0)
		{
			tempptr = currentptr;
			currentptr = currentptr->nextptr;
			delete tempptr;
		}
	}
	
}*/

bool EdgeList::isEmpty() const
{ return firstptr == 0; }

void EdgeList::insertnode(Edge* newptr) 
{
	if(isEmpty())
		firstptr = lastptr = newptr;
	else
	{
		lastptr->nextptr = newptr;
		lastptr = newptr;
	}
}

void EdgeList::print() const
{
	if(isEmpty())
	{
		cout<<"The list is empty.\n";
		return;
	}
	
	Edge* currentptr = firstptr;
	
	while(currentptr != 0)
	{
		cout<<currentptr<<endl;
		currentptr = currentptr->nextptr;
	}
}

void EdgeList::findPair()
{
	Vector v1a, v1b, v2a, v2b;
	Edge* currentptr1 = this->firstptr;
	Edge* currentptr2;

	while(currentptr1 != 0)
	{
		currentptr2 = currentptr1->nextptr;
		while(currentptr2 != 0)
		{
			v1a = currentptr1->verta->getVector();
			v1b = currentptr1->vertb->getVector();
			v2a = currentptr2->verta->getVector();
			v2b = currentptr2->vertb->getVector();
			
			if((myabs(v1a.x - v2b.x)<TOL) && (myabs(v1a.y - v2b.y)<TOL) && (myabs(v1a.z - v2b.z)<TOL) 
				&& (myabs(v1b.x - v2a.x)<TOL) && (myabs(v1b.y - v2a.y)<TOL) && (myabs(v1b.z - v2a.z)<TOL))
			{
				currentptr1->setPair(currentptr2);
				currentptr2->setPair(currentptr1);
			}
			else if((myabs(v1a.x - v2a.x)<TOL) && (myabs(v1a.y - v2a.y)<TOL) && (myabs(v1a.z - v2a.z)<TOL)
				&& (myabs(v1b.x - v2b.x)<TOL) && (myabs(v1b.y - v2b.y)<TOL) && (myabs(v1b.z - v2b.z)<TOL))
			{
 				currentptr2->verta->setVector(v2b);
				currentptr2->vertb->setVector(v2a);
				currentptr1->setPair(currentptr2);
				currentptr2->setPair(currentptr1);
			}
			currentptr2 = currentptr2->nextptr;
		}
		currentptr1 = currentptr1->nextptr;
	}
}
	
/******************************************************************************/

TriList::TriList():firstptr(0), lastptr(0) {}

TriList::~TriList()
{
	if(!isEmpty())
	{
		Triangle* currentptr = firstptr;
		Triangle* tempptr;
		
		while(currentptr != 0)
		{
			tempptr = currentptr;
			currentptr = currentptr->nextptr;
			delete tempptr;
		}
	}
	
}

bool TriList::isEmpty() const
{ return firstptr == 0; }

void TriList::insertnode(Triangle* newptr) 
{
	if(isEmpty())
	{
		firstptr = lastptr = newptr;
		firstptr->group_num = 1;
	}
	else
	{
		lastptr->nextptr = newptr;
		lastptr = newptr;
	}
}

void TriList::print() const
{
	if(isEmpty())
	{
		cout<<"The list is empty.\n";
		return;
	}
	
	Triangle* currentptr = firstptr;
	
	while(currentptr != 0)
	{
		cout<<currentptr<<endl;
		currentptr = currentptr->nextptr;
	}
}

void TriList::AvgNorm() const	// Calculating average normals without grouping
{
	Triangle* triptr = this->firstptr;
	
	float ctheta;
	Edge *edge, *edgestart, *pair;
	
	Triangle *face;
	Vector norm, norm1, avg_norm;
	Vector sum;
	Vertex* vert;

	while(triptr != 0)
	{
		for(int i = 0; i<3; i++)
		{
			int count = 1;
			vert = triptr->vert[i];
			edge = vert->getEdge();
			edgestart = edge;
			
			norm1 = triptr->facet_norm;
			sum = norm1;

			do
			{
				pair = edge->getPair();
				if(pair == 0)
				{
					cout<<"Error in STL file!"<<endl;
					exit(0);
				}
				edge = pair->getNext();
				face = edge->getFace();
				norm = face->facet_norm;
				ctheta = norm & norm1;

				//if(ctheta > (1 - TOL))
				if(ctheta > 0.90)
				{
					sum = sum + norm;
					count++;
				}
			}
			while(edge != edgestart);
	
			avg_norm = sum.scale(1.0/(float) count);
			vert->setAvgNorm(avg_norm);
			
			//cout<<avg_norm.x<<" "<<avg_norm.y<<" "<<avg_norm.z<<" "<<avg_norm.mod()<<endl;
		}
		triptr = triptr->nextptr;
	}
}

// Gives better smoothing than non- grouping method
void TriList::computeGroup(Triangle* currentTri)	// Computing group #s
{
	int flag = 0;
	Triangle* adjface;
	Edge* edgeStart = currentTri->edge;
	Edge* edgeCount = currentTri->edge;
	Vector currentNorm = currentTri->facet_norm;
	Vector v1, v2;
	Vector adjNorm;
	Edge* pair;
	float threshold;

	if(currentTri->bits == 3)
		return;
	for(int i=0; i<3; i++)
	{
		pair = edgeCount->getPair();
		if(pair == 0)
		{
			cout<<"Error in STL file!"<<endl;
			exit(0);
		}
		adjface = pair->getFace();
		adjface->bits++;
		edgeCount = edgeCount->getNext();
	
		if(adjface->group_num > 0)
		{
			//flag = 1;
			continue;
		}
		else
		{
			adjNorm = adjface->facet_norm;
			// Normalizing facet normals
			v1 = ~adjNorm;		
			v2 = ~currentNorm;
			threshold = v1&v2;
	
			//if(threshold < (1 - TOL))
			if(threshold < 0.95)	
				adjface->group_num = currentTri->group_num + 1;
			else 
				adjface->group_num = currentTri->group_num;			
		}
	
		computeGroup(adjface);
	}	// end for
}

