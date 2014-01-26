#include "group.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

extern int objcount;
extern genList<Object> objectlist;
extern int nx, ny;
extern float sx, sy;
extern Vector eye, aim;
//extern Vector light[5];
extern Light light[5];
extern Color ibackground;

const int DEPTH = 4;
const int LIGHT = 1;

int Partition(ObjCopy* , int, int);

void QuickSort(ObjCopy* A, int p, int r)
{
	int q;
	if(p < r)
	{
		q = Partition(A, p, r);
		QuickSort(A, p, q);
		QuickSort(A, q+1, r);
	}
}

int Partition(ObjCopy* A, int p, int r)
{
	ObjCopy temp;
	float x = A[p].tvalue;
	int i = p-1;
	int j = r+1;

	while(1)
	{
		do{
			j--;
		}
		while(A[j].tvalue > x);
		
		do{
			i++;
		}
		while(A[i].tvalue < x);

		if(i < j)
		{
			temp = A[i];
			A[i] = A[j];
			A[j] = temp;
		}
		else
			return j;
	}
}

// Check if return value = -1 or < 0.001 => invalid intersection
// face normal (normalized), triangle vertex point
float Ray::intersect_plane(Vector faceNormal, Vector point)
{
	Vector vsub;
	float dot1, dot2;
	float t;

	vsub = point - this->base;
	dot1 = faceNormal & vsub;
	dot2 = faceNormal & this->dir;
	if(myabs(dot2 - 0.0) <= TOL)
		return -1;
	else
	{
		t = dot1/dot2;
		return t;
	}
	
}

// returns 0 if point is outside polygon, 1 otherwise
// v1, v2, v3 are the vertices of each triangle
// Vector Q = ray.base + ray.dir.scale(t);
int pip(Vector v1, Vector v2, Vector v3, Vector Q, 
			 Vector faceNormal /* normalized */)	// function for Point in Polygon (PIP) test
{
	Vector vsub, vcross, edge;
	float vdot1, vdot2, vdot3;

	edge = v2 - v1;
	vsub = Q - v1;
	vcross = edge * vsub;
	vdot1 = vcross & faceNormal;

	edge = v3 - v2;
	vsub = Q - v2;
	vcross = edge * vsub;
	vdot2 = vcross & faceNormal;

	edge = v1 - v3;
	vsub = Q - v3;
	vcross = edge * vsub;
	vdot3 = vcross & faceNormal;

	if(((vdot1 >= 0)&&(vdot2 >= 0)&&(vdot3 >= 0)) || ((vdot1 <= 0)&&(vdot2 <= 0)&&(vdot3 <= 0)))
		return 1; 
	else
		return 0;
}

// n1, n2 and n3 are the three avg normals at the vertices of each triangle
Vector Triangle::Blend_norm(Vector Q)
{
	Vector v1, v2, v3;
	Vector n1, n2, n3;
	Vector vsub1, vsub2;
	Vector a1, a2, a3;
	Vector blend;
	float area1, area2, area3, sum;
	float t1, t2, t3;

	n1 = this->vert[0]->getAvgNorm();
	n2 = this->vert[1]->getAvgNorm();
	n3 = this->vert[2]->getAvgNorm();
	v1 = this->vert[0]->getVector();
	v2 = this->vert[1]->getVector();
	v3 = this->vert[2]->getVector();

	vsub1 = Q - v2;
	vsub2 = v3 - v2;
	a1 = vsub1*vsub2;

	vsub1 = Q - v3;
	vsub2 = v1 - v3;
	a2 = vsub1*vsub2;

	vsub1 = Q - v1;
	vsub2 = v2 - v1;
	a3 = vsub1*vsub2;

	area1 = a1.mod();
	area2 = a2.mod();
	area3 = a3.mod();

	sum = area1+area2+area3;
	t1 = area1/sum;
	t2 = area2/sum;
	t3 = area3/sum;

	blend = n1.scale(t1) + n2.scale(t2) + n3.scale(t3);
	return ~blend;
}


float Object::bound_box(Ray ray) const
{
	float t, thit = HUGEREAL;
	float xmin, xmax, ymin, ymax, zmin, zmax;
	Vector Q;
	Vector i(1,0,0), j(0,1,0), k(0,0,1);
	Vector base = ray.getBase();
	Vector dir = ray.getDir();
	
	xmin = this->vmin.x;
	xmax = this->vmax.x;
	ymin = this->vmin.y;
	ymax = this->vmax.y;
	zmin = this->vmin.z;
	zmax = this->vmax.z;

	t = ray.intersect_plane(i, this->vmax);
	if(t > 0.001)
	{
		Q = base + dir.scale(t);

		if( (ymin <= Q.y) && (Q.y <= ymax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(i.scale(-1), this->vmin);
	if(t > 0.001)
	{
		Q = base + dir.scale(t);
		
		if( (ymin <= Q.y) && (Q.y <= ymax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(j, this->vmax);
	if(t > 0.001)
	{
		Q = base + dir.scale(t);
				
		if( (xmin <= Q.x) && (Q.x <= xmax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(j.scale(-1), this->vmin);
	if(t > 0.001)
	{
		Q = base + dir.scale(t);
		
		if( (xmin <= Q.x) && (Q.x <= xmax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(k, this->vmax);
	if(t > 0.001)
	{
		Q = base + dir.scale(t);
		
		if( (xmin <= Q.x) && (Q.x <= xmax) && (ymin <= Q.y) && (Q.y <= ymax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(k.scale(-1), this->vmin);
	if(t > 0.001)
	{
		Q = base + dir.scale(t);
		
		if( (xmin <= Q.x) && (Q.x <= xmax) && (ymin <= Q.y) && (Q.y <= ymax) )
		{
			thit = t;
			return thit;
		}
	}
	return thit;
}
		
float Object::Shadow_box(Ray ray) const	// Bounding box function for shadows
{
	float t, thit = 1.0;
	float xmin, xmax, ymin, ymax, zmin, zmax;
	Vector Q;
	Vector i(1,0,0), j(0,1,0), k(0,0,1);
	Vector base = ray.getBase();
	Vector dir = ray.getDir();
	
	xmin = this->vmin.x;
	xmax = this->vmax.x;
	ymin = this->vmin.y;
	ymax = this->vmax.y;
	zmin = this->vmin.z;
	zmax = this->vmax.z;

	t = ray.intersect_plane(i, this->vmax);
	if( (t > 0.0) && (t < 1.0) )
	{
		Q = base + dir.scale(t);

		if( (ymin <= Q.y) && (Q.y <= ymax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(i.scale(-1), this->vmin);
	if( (t > 0.0) && (t < 1.0) )
	{
		Q = base + dir.scale(t);
		
		if( (ymin <= Q.y) && (Q.y <= ymax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(j, this->vmax);
	if( (t > 0.0) && (t < 1.0) )
	{
		Q = base + dir.scale(t);
				
		if( (xmin <= Q.x) && (Q.x <= xmax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(j.scale(-1), this->vmin);
	if( (t > 0.0) && (t < 1.0) )
	{
		Q = base + dir.scale(t);
		
		if( (xmin <= Q.x) && (Q.x <= xmax) && (zmin <= Q.z) && (Q.z <= zmax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(k, this->vmax);
	if( (t > 0.0) && (t < 1.0) )
	{
		Q = base + dir.scale(t);
		
		if( (xmin <= Q.x) && (Q.x <= xmax) && (ymin <= Q.y) && (Q.y <= ymax) )
		{
			thit = t;
			return thit;
		}
	}
	t = ray.intersect_plane(k.scale(-1), this->vmin);
	if( (t > 0.0) && (t < 1.0) )
	{
		Q = base + dir.scale(t);
		
		if( (xmin <= Q.x) && (Q.x <= xmax) && (ymin <= Q.y) && (Q.y <= ymax) )
		{
			thit = t;
			return thit;
		}
	}
	return thit;
}

// Reflect_ray() and Transmit_ray() return unit directions

// normalized, normalized (blended normal)
Ray Ray::Reflect_ray() 
{
	Ray ray;
	Vector norm = this->blend_norm;
	float dot;
	/* Normal Flipping */
	if( (this->dir & norm) > 0.0 )
		norm = norm.scale(-1);
	
	dot = this->dir & norm;
	ray.dir = ~(norm.scale(-2*dot) + this->dir);
	ray.base = this->hitpoint;
	return ray;
}

// normalized blended normal
Ray Ray::Transmit_ray()
{
	Ray ray;
	float index;
	float dot, factor, ratio;
	Vector v1, v2, v3, norm;
	
	index = this->object->getIndex();
	norm = this->blend_norm;

	/* Normal Flipping */
	if( (this->dir & norm) < 0.0)
	{
		ratio = 1.0/index;
	}
	else
	{
		ratio = index;
		norm = norm.scale(-1);
	}
		

	dot = this->dir & norm;
	v1 = norm.scale(-1*dot) + this->dir;
	v2 = v1.scale(ratio);

	/* Accounting for total internal reflection */
	if( (1 - (ratio*ratio* (1 - dot*dot)) ) < 0.0)
		ray = this->Reflect_ray();
	else
	{
		factor = pow((1 - (ratio*ratio* (1 - dot*dot))), 0.5);
		v3 = norm.scale(factor);
		ray.base = this->hitpoint;
		ray.dir = ~(v2 - v3);
	}

	factor = pow((1 - (ratio*ratio* (1 - dot*dot))), 0.5);
	v3 = norm.scale(factor);
	ray.base = this->hitpoint;
	ray.dir = ~(v2 - v3);
	return ray;
}

// returns intensity
// eye is global property
Color Ray::Shadow_ray()
{
	Color ishadow;
	Object* current_obj;
	Triangle *current_tri, *lastTri;
	Vertex** Vert;
	Vector vmin, vmax;
	Ray shadow;
	Vector Q;
	int p;
	float t, thit = 1.0;
	
	//float sumd = 0.0, sums = 0.0;
	float sumd_r = 0.0, sumd_g = 0.0, sumd_b = 0.0, sums_r = 0.0, sums_g = 0.0, sums_b = 0.0;
	float ctheta, cphi, cnphi;
	Color Kd, iamb;
	float Kd_r, Kd_g, Kd_b;
	float Ks;
	float iamb_r, iamb_g, iamb_b, n;
	Vector norm = this->blend_norm;	// blended normal
	
	if( (this->dir & norm) > 0.0 )
		norm = norm.scale(-1);
	
	//Vector norm = this->triangle->getnormal();	// flat shading
	Vector h;
	Vector lightray, eyeray;	
	
	for(int i=0; i<LIGHT; i++)
	{
		//lightray = light[i] - this->hitpoint;
		lightray = light[i].v - this->hitpoint;
		
		eyeray = eye - this->hitpoint;
		ctheta = norm & (~lightray);
		if(ctheta < 0)
			ctheta = 0;
		h = ~(eyeray + lightray);
		cphi = norm & h;
				
		n = this->object->getN();
		cnphi = pow(cphi, n);
		if(cnphi < 0)
			cnphi = 0;

		shadow.base = this->hitpoint;
		shadow.dir = lightray;

		current_obj = objectlist.firstptr;
		while(current_obj != 0)
		{
			t = current_obj->Shadow_box(shadow);
			if(t == 1)
			{
				current_obj = current_obj->nextptr;
				continue;
			}
			else
			{
				current_tri = current_obj->getFirstTri();
				lastTri = current_obj->getLastTri();
				while(current_tri != lastTri->nextptr)
				{
					Vert = current_tri->getVerts();
					t = shadow.intersect_plane(current_tri->getnormal(), Vert[0]->getVector());
					if((t > TOL) && (t < (1.0 - TOL)))
					{
						Q = shadow.base + shadow.dir.scale(t);
						p = pip(Vert[0]->getVector(), Vert[1]->getVector(), 
								   Vert[2]->getVector(), Q, current_tri->getnormal());
						if(p == 1)
						{
							ctheta = 0;
							cnphi = 0;
							break;
						}
					}
					current_tri = current_tri->nextptr;
				}
			}
			current_obj = current_obj->nextptr;
		}
		/*sumd = sumd + ctheta;
		sums = sums + cnphi;*/
		sumd_r = sumd_r + light[i].col.red*ctheta;
		sumd_g = sumd_g + light[i].col.green*ctheta;
		sumd_b = sumd_b + light[i].col.blue*ctheta;
		
		sums_r = sums_r + light[i].col.red*cnphi;
		sums_g = sums_g + light[i].col.green*cnphi;
		sums_b = sums_b + light[i].col.blue*cnphi;
	}
	
	Kd = this->object->getKdiff();
	Kd_r = Kd.red;
	Kd_g = Kd.green;
	Kd_b = Kd.blue;
	Ks = this->object->getKspec(); 
	
	iamb = this->object->getIamb();
	iamb_r = iamb.red*light[i].col.red;
	iamb_g = iamb.green*light[i].col.green;
	iamb_b = iamb.blue*light[i].col.blue;
	
	/*ishadow.red = iamb_r + Kd_r*sumd + Ks*sums;	
	ishadow.green = iamb_g + Kd_g*sumd + Ks*sums;	
	ishadow.blue = iamb_b + Kd_b*sumd + Ks*sums;	*/

	ishadow.red = iamb_r + Kd_r*sumd_r + Ks*sums_r;	
	ishadow.green = iamb_g + Kd_g*sumd_g + Ks*sums_g;	
	ishadow.blue = iamb_b + Kd_b*sumd_b + Ks*sums_b;	
	return ishadow;
}

Color Intensity(Ray ray)	// Ray Tracing Recursion
{
	/*Color WhiteLight;
	WhiteLight.red = 1.0;
	WhiteLight.blue = 1.0;
	WhiteLight.green = 1.0;*/

	int depth = ray.getDepth();
	Vector lightray;
	int p;
	Color Ish, Iref, Itrans, I;
	Color Kr, Kt;
	Ray ray1, ray2;
	float Kr_r, Kr_g, Kr_b;
	float Kt_r, Kt_g, Kt_b;
	Object* obj;

	if(depth > DEPTH)
			return ibackground;	// Background intensity (global)
	else	// 1
	{
		p = ray.Compute_intersect();
		if(p == 0)
		{
			/*for(int i=0; i<LIGHT; i++)
			{
				lightray = ~(light[i] - ray.getBase());
				
				if( (1.0 - (ray.getDir() & lightray)) < TOL)
					return WhiteLight;	 //return <light[i].intensity>
				else
					return ibackground;
			}*/
			return ibackground;
		}
		else	// 2
		{
			obj = ray.getObject();
	
			Ish = ray.Shadow_ray();
			if(Ish.red > 1.0)
				Ish.red = 1.0;
			if(Ish.blue > 1.0)
				Ish.blue = 1.0;
			if(Ish.green > 1.0)
				Ish.green = 1.0;
						
			/*Kr = obj->getKref();
			Kr_r = Kr.red;
			Kr_g = Kr.green;
			Kr_b = Kr.blue;

			if( (Kr_r == 0.0) && (Kr_g == 0.0) && (Kr_b == 0.0) )
			{
				Iref.red = 0.0;
				Iref.blue = 0.0;
				Iref.green = 0.0;
			}
			else
			{
				ray1 = ray.Reflect_ray();
				ray1.setDepth(depth + 1);
				Iref = Intensity(ray1);	// Recursive call
			
				if(Iref.red > 1.0)
					Iref.red = 1.0;
				if(Iref.blue > 1.0)
					Iref.blue = 1.0;
				if(Iref.green > 1.0)
					Iref.green = 1.0;
			}

			Kt = obj->getKtrans();
			Kt_r = Kt.red;
			Kt_g = Kt.green;
			Kt_b = Kt.blue;

			if( (Kt_r == 0.0) && (Kt_g == 0.0) && (Kt_b == 0.0) )
			{
				Itrans.red = 0.0;
				Itrans.blue = 0.0;
				Itrans.green = 0.0;
			}
			else
			{
				ray2 = ray.Transmit_ray();
				ray2.setDepth(depth + 1);

				Itrans = Intensity(ray2);	// Recursive call
				
				if(Itrans.red > 1.0)
					Itrans.red = 1.0;
				if(Itrans.blue > 1.0)
					Itrans.blue = 1.0;
				if(Itrans.green > 1.0)
					Itrans.green = 1.0;
			}
			
			I.red = Ish.red + Kr_r*Iref.red + Kt_r*Itrans.red;
			I.green = Ish.green + Kr_g*Iref.green + Kt_g*Itrans.green;
			I.blue = Ish.blue + Kr_b*Iref.blue + Kt_b*Itrans.blue;
			return I;*/
			return Ish;

		}	// end else 2
	}	// end else 1
}

  
// int nx, int ny, int sx, int sy, Vector eye, Vector aim are all global
void Ray_trace()
{
	Ray parent;
	Vector ze = ~(eye - aim);
	Vector up(0,1,0);
	Vector xebar, xe, ye;
	float dx, dy;
	Color Out;
	int i, j;

	cout<<"# of objects = "<<objcount<<endl;
	ofstream out_file("test.ppm", ios::out);
	out_file<<"P3\n";
	out_file<<nx<<" "<<ny<<endl;
	out_file<<255<<endl;

	xebar = up*ze;
	if(myabs(xebar.mod() - 0.0) < TOL)
	{
		if(ze.y > 0)
		{
			up.x = 0;
			up.y = 0;
			up.z = -1;
		}
		else
		{
			up.x = 0;
			up.y = 0;
			up.z = 1;
		}
		xebar = up*ze;
	}
	xe = ~xebar;
	ye = ze*xe;

	Vector s, sfirst;
	Vector rd;

	sfirst = aim - xe.scale(sx/2.0) + ye.scale(sy/2.0);
	s = sfirst;
	
	dx = (float)sx/(nx-1);
	dy = (float)sy/(ny-1);

	for(i=0; i<ny; i++)
	{
		for(j=0; j<nx; j++)
		{
			rd = ~(s - eye);
			parent.setBase(eye);
			parent.setDir(rd);
			parent.setDepth(0);
			Out = Intensity(parent);	// parent ray through each pixel

			int x = (int)255*Out.red;
			int y = (int)255*Out.green;
			int z = (int)255*Out.blue;
			if(x > 255)
				x = 255;
			if(y > 255)
				y = 255;
			if(z > 255)
				z = 255;
				
			//cout<<x<<" "<<y<<" "<<z<<" ";
			out_file<<x<<" "<<y<<" "<<z<<" ";
			s = s + xe.scale(dx);
		}
		out_file<<"\n";
		s = sfirst - ye.scale(dy*(i+1));
	}
	out_file.close();
}

int Ray::Compute_intersect()
{
	Color In;
	In.red = 1.0;
	In.green = 0.0;
	In.blue = 0.0;

	ObjCopy* copy = new ObjCopy[objcount];
	float t, thit = HUGEREAL, thitTri = HUGEREAL;
	int p;
	int count = 0;
	Vector v1, v2, v3, n1, n2, n3;
	Vector blend, Q;
	Object *curr_obj, *objecthit;
	Triangle *curr_tri, *lastTri;
	Triangle* trihit = 0;
	Vertex** Vert;

	Vector base = this->base;
	Vector dir = this->dir;

	curr_obj = objectlist.firstptr;
	while(curr_obj != 0)
	{
		t = curr_obj->bound_box(*this);
		copy[count].obj = curr_obj;
		copy[count].tvalue = t;
		count++;
		curr_obj = curr_obj->nextptr;
	}
	
	QuickSort(copy, 0, objcount-1);
	
	for(count = 0; count < objcount; count++)
	{
		thit = copy[count].tvalue;
		objecthit = copy[count].obj;
	
		if(thit == HUGEREAL)
			return 0;
		else
		{
			curr_tri = objecthit->getFirstTri();
			lastTri = objecthit->getLastTri();
			
			while(curr_tri != lastTri->nextptr)
			{
				Vert = curr_tri->getVerts();
				v1 = Vert[0]->getVector();
				v2 = Vert[1]->getVector();
				v3 = Vert[2]->getVector();
	
				t = this->intersect_plane(curr_tri->getnormal(), v1);
						
				if(t > 0.001) 
				{
					Q = base + dir.scale(t);
					p = pip(v1, v2, v3, Q, curr_tri->getnormal());
					
					if(p == 1) 
					{
						if(t < thitTri)
						{
							thitTri = t;
							trihit = curr_tri;						
						}
					}	
				}
				curr_tri = curr_tri->nextptr;
				
			}
			if(trihit != 0)
			{
				Q = base + dir.scale(thitTri);
				this->hitpoint = Q;
				blend = trihit->Blend_norm(Q);
				this->triangle = trihit;
				this->blend_norm = blend;
				this->object = objecthit;
				//cout<<blend.x<<" "<<blend.y<<" "<<blend.z<<" "<<blend.mod()<<endl;
				return 1;
	
			}
			
		}
	}	
	return 0;
}

// No Bounding Box
/*Color Ray::Compute_intersect()	// Compare time using bounding box and without it
{
	Color In;
	In.red = 1.0;
	In.green = 0.0;
	In.blue = 0.0;
	
	float t, thit = HUGEREAL, thitTri = HUGEREAL;
	int p;
	Vector v1, v2, v3, n1, n2, n3;
	Vector blend, Q;
	Object *curr_obj, *objecthit;
	Triangle *curr_tri, *lastTri;
	Triangle* trihit = 0;
	Vertex** Vert;

	Vector base = this->base;
	Vector dir = this->dir;

	curr_obj = objectlist.firstptr;

	while(curr_obj != 0)	
	{	
		curr_tri = curr_obj->getFirstTri();	
		lastTri = curr_obj->getLastTri();	

		while(curr_tri != lastTri->nextptr)
		{
			Vert = curr_tri->getVerts();
			v1 = Vert[0]->getVector();
			v2 = Vert[1]->getVector();
			v3 = Vert[2]->getVector();

			t = this->intersect_plane(curr_tri->getnormal(), v1);
			//cout<<t<<endl;
			if(t >0.001)	
			{
				//cout<<t<<endl;
				Q = base + dir.scale(t);
				p = pip(v1, v2, v3, Q, curr_tri->getnormal());
				if(p == 1)
				{
					if(t < thit)
					{
						thit = t;
						trihit = curr_tri;
						objecthit = curr_obj;
					}
				}
			}

			curr_tri = curr_tri->nextptr;
		}
		curr_obj = curr_obj->nextptr;
		
	}
	//cout<<thit<<" "<<objecthit<<endl;
	if(thit == HUGEREAL)
		return ibackground;
	else
	{
		Q = base + dir.scale(thit);
		this->hitpoint = Q;
		blend = trihit->Blend_norm(Q);
		this->depth++;
		this->blend_norm = blend;
		cout<<blend.x<<" "<<blend.y<<" "<<blend.z<<" "<<blend.mod()<<endl;
		this->object = curr_obj;

		return In;

	}
}*/
