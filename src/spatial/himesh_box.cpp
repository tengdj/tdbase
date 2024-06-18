/*
 * himesh_box.cpp
 *
 *  Created on: Aug 30, 2023
 *      Author: teng
 */

#include "aab.h"
#include <unistd.h>
#include <cstring>
#include "util.h"

using namespace std;

namespace tdbase{

/*
 * functions for the aab class
 *
 * */

aab::aab(){
	reset();
}

aab::aab(const aab &b){
	for(int i=0;i<3;i++){
		low[i] = b.low[i];
		high[i] = b.high[i];
	}
}
aab::aab(float min_x, float min_y, float min_z,
		float max_x, float max_y, float max_z){
	low[0] = min_x;
	low[1] = min_y;
	low[2] = min_z;
	high[0] = max_x;
	high[1] = max_y;
	high[2] = max_z;
}
void aab::set_box(float l0, float l1, float l2, float h0, float h1, float h2){
	low[0] = l0;
	low[1] = l1;
	low[2] = l2;
	high[0] = h0;
	high[1] = h1;
	high[2] = h2;
}
void aab::reset(){
	for(int i=0;i<3;i++){
		low[i] = DBL_MAX;
		high[i] = -DBL_MAX;
	}
}

void aab::update(float x, float y, float z){
	if(low[0]>x){
		low[0] = x;
	}
	if(high[0]<x){
		high[0] = x;
	}
	if(low[1]>y){
		low[1] = y;
	}
	if(high[1]<y){
		high[1] = y;
	}
	if(low[2]>z){
		low[2] = z;
	}
	if(high[2]<z){
		high[2] = z;
	}
}

void aab::update(const aab &p){
	for(int i=0;i<3;i++){
		if(low[i]>p.low[i]){
			low[i] = p.low[i];
		}
		if(high[i]<p.high[i]){
			high[i] = p.high[i];
		}
	}
}
void aab::set_box(const aab &b){
	for(int i=0;i<3;i++){
		low[i] = b.low[i];
		high[i] = b.high[i];
	}
}


bool aab::intersect(aab &object){
	return !(object.low[0] >= high[0] || object.high[0] <= low[0] ||
			 object.low[1] >= high[1] || object.high[1] <= low[1] ||
			 object.low[2] >= high[2] || object.high[2] <= low[2]);
}
bool aab::contains(aab *object){
	for(int i=0;i<3;i++){
		if(object->low[i]<low[i]){
			return false;
		}
		if(object->high[i]>high[i]){
			return false;
		}
	}
	return true;
}

bool aab::contains(float *point){
	for(int i=0;i<3;i++){
		if(point[i]<low[i]){
			return false;
		}
		if(point[i]>high[i]){
			return false;
		}
	}
	return true;
}

void aab::print(){
	cout<<"(";
	cout<<low[0]<<",";
	cout<<low[1]<<",";
	cout<<low[2]<<")";
	cout<<" -> (";
	cout<<high[0]<<",";
	cout<<high[1]<<",";
	cout<<high[2]<<")";
	cout<<endl;
}

float aab::diagonal_length(){
	float dl = 0;
	for(int i=0;i<3;i++){
		dl += (high[i]-low[i])*(high[i]-low[i]);
	}
	return sqrt(dl);
}
float aab::volume(){
	return (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
}

// get the possible minimum and maximum distance of
// objects with their aabs
range aab::distance(const aab &b){
	range ret;
	ret.maxdist = 0;
	ret.mindist = 0;
	float tmp1 = 0;
	float tmp2 = 0;
	for(int i=0;i<3;i++){
		tmp1 = low[i]-b.high[i];
		tmp2 = high[i]-b.low[i];
		ret.maxdist += (tmp1+tmp2)*(tmp1+tmp2)/4;
		if(tmp2<0){
			ret.mindist += tmp2*tmp2;
		}else if(tmp1>0){
			ret.mindist += tmp1*tmp1;
		}
	}
	ret.mindist = sqrt(ret.mindist);
	ret.maxdist = sqrt(ret.maxdist);
	return ret;
}


/*
 *
 * Voxel functions
 *
 * */

Voxel::~Voxel(){
	clear();
}

// clear the buffer
void Voxel::clear(){
	if(owned && triangles){
		delete []triangles;
		triangles = NULL;
	}
	if(owned && hausdorff){
		delete []hausdorff;
		hausdorff = NULL;
	}
	capacity = 0;
	owned = false;
}

// create buffer if needed
void Voxel::reserve(size_t size){
	num_triangles = 0;
	if(size<=capacity){
		return;
	}
	clear();
	capacity = size;
	triangles = new float[9*size];
	hausdorff = new float[2*size];
	owned = true;
}

void Voxel::insert(const float *t, const float *h){
	// enough space must be researved
	assert(num_triangles < capacity);
	memcpy(triangles+num_triangles*9, t, 9*sizeof(float));
	memcpy(hausdorff+num_triangles*2, h, 2*sizeof(float));
	num_triangles++;
}

void Voxel::print(){
	aab::print();
	if(num_triangles == 0){
		printf("voxel is empty\n");
		return;
	}
	for(int i=0;i<num_triangles;i++){
		float *tri = triangles + i*9;
		printf("(");
		for(int j=0;j<9;j++){
			if(j==3 || j==6){
				printf(",");
			}
			printf("%.2f ", tri[j]);
		}
		printf(")|(%.3f, %.3f)\n",*(hausdorff+i*2),*(hausdorff+i*2+1));
	}
}

// copy to the buffer
void Voxel::batch_load(const float *t, const float *h, size_t s){
	reserve(s);
	memcpy(triangles, t, s*9*sizeof(float));
	memcpy(hausdorff, h, s*2*sizeof(float));
	num_triangles = s;
//	print();
//	log("");
}

// link the buffer to external spaces
void Voxel::external_load(float *t, float *h, size_t s){
	clear();
	// should not be freed in the deconstructor
	triangles = t;
	hausdorff = h;
	num_triangles = s;
}

float Voxel::getHausdorffDistance(int offset){
	assert(hausdorff);
	assert(offset>=0 && offset < num_triangles);
	return *(hausdorff+offset*2+1);
}

float Voxel::getProxyHausdorffDistance(int offset){
	assert(hausdorff);
	assert(offset>=0 && offset < num_triangles);
	return *(hausdorff+offset*2);
}


}




