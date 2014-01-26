#include "group.h"
#include<iostream>
#include<fstream>

using namespace std;

const int SIZE = 255;

extern genList<Object> objectlist;
Object* item;
int count;
int objcount = 0;
int nx, ny;
float sx, sy;
float n;
float index;
Vector eye;
Vector aim;
//Vector light[5];
Light light[5];
Color ibackground;

void Readprop(const char* filename)
{
	count = 0;
	string keyword;
	
	ifstream input_file(filename, ios::in);
	if(!input_file.is_open())
	{
		cout <<"Error: Cannot input filename: "<< filename << endl;
		exit(0);
	}

	item = objectlist.firstptr;
	while(input_file.peek() != EOF)
	{
		int j = 0;
		char c;
		keyword = "";
		char line[SIZE];

		while(isalnum(c = input_file.get()))
		{
			line[j] = c;
			j++;
		}
		line[j] = NULL;
		keyword = line;
		//cout<<keyword<<endl;

		if((keyword == "NX") || (keyword == "nx"))
		{
			input_file >> nx;
			input_file.getline(line, SIZE, '\n');
			cout<<nx<<endl;
		}
		
		if((keyword == "NY") || (keyword == "ny"))
		{
			input_file >> ny;
			input_file.getline(line, SIZE, '\n');
			cout<<ny<<endl;
		}
		
		if((keyword == "SX") || (keyword == "sx"))
		{
			input_file >> sx;
			input_file.getline(line, SIZE, '\n');
			cout<<sx<<endl;
		}
		
		if((keyword == "SY") || (keyword == "sy"))
		{
			input_file >> sy;
			input_file.getline(line, SIZE, '\n');
			cout<<sy<<endl;
		}
		

		if((keyword == "EYEPOINT") || (keyword == "eyepoint"))
		{
			input_file >> eye.x >> eye.y >> eye.z;
			input_file.getline(line, SIZE, '\n');
			cout<<eye.x<<" "<<eye.y<<" "<<eye.z<<endl;
		}

		if((keyword == "AIMPOINT") || (keyword == "aimpoint"))
		{
			input_file >> aim.x >> aim.y >> aim.z;
			input_file.getline(line, SIZE, '\n');
			cout<<aim.x<<" "<<aim.y<<" "<<aim.z<<endl;
		}

		if((keyword == "LIGHTSOURCE") || (keyword == "lightsource"))
		{
			//input_file >> light[count].x >> light[count].y >> light[count].z;
			input_file >> light[count].v.x >> light[count].v.y >> light[count].v.z>>light[count].col.red >> light[count].col.green>>light[count].col.blue;
			input_file.getline(line, SIZE, '\n');
			//cout<<light[count].x<<" "<<light[count].y<<" "<<light[count].z<<endl;
			count++;
		}
				

		if((keyword == "BACKGROUND") || (keyword == "background"))
		{
			input_file>> ibackground.red >> ibackground.green >> ibackground.blue;
			input_file.getline(line, SIZE, '\n');
			cout<<ibackground.red<<" "<<ibackground.green<<" "<<ibackground.blue<<endl;
		}
			
		if((keyword == "PROPERTY") || (keyword == "property"))
		{
			//input_file.getline(line, SIZE, '\n');
			objcount++;
			cout<<item<<endl;
		}
			
		if((keyword == "IAMBIENT") || (keyword == "iambient"))
		{
			Color iamb, Iamb;
			input_file >> iamb.red >> iamb.green >> iamb.blue;
			cout<<iamb.red<<" "<<iamb.green<<" "<<iamb.blue<<endl;
			item->setIamb(iamb);
			input_file.getline(line, SIZE, '\n');
			Iamb = item->getIamb();
			cout<<Iamb.red<<" "<<Iamb.green<<" "<<Iamb.blue<<endl;
		}

		if((keyword == "KSPECULAR") || (keyword == "kspecular"))
		{
			float Ks, K;
			input_file >> Ks;
			item->setKspec(Ks);
			input_file.getline(line, SIZE, '\n');
			K = item->getKspec();
			cout<<K<<endl;
		}

		if((keyword == "KDIFFUSED") || (keyword == "kdiffused"))
		{
			Color Kd, K;
			input_file >> Kd.red >> Kd.green >> Kd.blue;
			item->setKdiff(Kd);
			input_file.getline(line, SIZE, '\n');
			K = item->getKdiff();
			cout<<K.red<<" "<<K.green<<" "<<K.blue<<endl;
		}

		if((keyword == "KREFLECT") || (keyword == "kreflect"))
		{
			Color Kr, K;
			input_file >> Kr.red >> Kr.green >> Kr.blue;
			item->setKref(Kr);
			input_file.getline(line, SIZE, '\n');
			K = item->getKref();
			cout<<K.red<<" "<<K.green<<" "<<K.blue<<endl;
		}

		if((keyword == "KTRANSMIT") || (keyword == "ktransmit"))
		{
			Color Kt, K;
			input_file >> Kt.red >> Kt.green >> Kt.blue;
			item->setKtrans(Kt);
			input_file.getline(line, SIZE, '\n');
			K = item->getKtrans();
			cout<<K.red<<" "<<K.green<<" "<<K.blue<<endl;
		}

		if((keyword == "INDEX") || (keyword == "index"))
		{
			input_file >> index;
			item->setIndex(index);
			input_file.getline(line, SIZE, '\n');
			cout<<item->getIndex()<<endl;
		}

		if((keyword == "N") || (keyword == "n"))
		{
			input_file >> n;
			item->setN(n);
			input_file.getline(line, SIZE, '\n');
			cout<<item->getN()<<endl;
		}

		if((keyword == "ENDPROPERTY") || (keyword == "endproperty"))
		{
			if(item->nextptr != 0)
				item = item->nextptr;
		}
	}
}
