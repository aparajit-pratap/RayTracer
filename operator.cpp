#include "group.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

Vector Vector::scale(float t)
{
	Vector product;
	product.x = this->x*t;
	product.y = this->y*t;
	product.z = this->z*t;
	return product;
}

Vector Vector::operator-(const Vector& v)
{

	Vector subv;
	subv.x = this->x - v.x;
	subv.y = this->y - v.y;
	subv.z = this->z - v.z;
	return subv;
}

Vector Vector::operator+(const Vector& v)
{

	Vector addv;
	addv.x = this->x + v.x;
	addv.y = this->y + v.y;
	addv.z = this->z + v.z;
	return addv;
}

Vector Vector::operator*(const Vector& v)
{
	Vector crossv;
	crossv.x = this->y*v.z - this->z*v.y;
	crossv.y = this->z*v.x - this->x*v.z;
	crossv.z = this->x*v.y - this->y*v.x;
	return crossv;
}

Vector Vector::operator~()
{
	Vector unitv;
	double den = sqrt((this->x)*(this->x) + (this->y)*(this->y) + (this->z)*(this->z));
	unitv.x = this->x/(float)den;
	unitv.y = this->y/(float)den;
	unitv.z = this->z/(float)den;
	return unitv;
}

float Vector::operator&(const Vector& v)
{
	float dotv;
	dotv = (this->x)*v.x + (this->y)*v.y + (this->z)*v.z;
	return dotv;
}
	
float Vector::mod()
{
	float modv;
	modv = (float)sqrt(this->x*this->x + this->y*this->y + this->z*this->z);
	return modv;
}

ostream& operator<<(ostream& out, const Vector& v)
{
	out <<v.x<<"i + "<<v.y<<"j + "<<v.z<<"k"<<endl;
	return out;
}

