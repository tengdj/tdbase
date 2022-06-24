/*
 * aab.h
 *
 *  Created on: Nov 14, 2019
 *      Author: teng
 *
 *  define of Axis-Aligned Box, which is a basic unit
 *  for describing objects
 */

#ifndef HISPEED_AAB_H_
#define HISPEED_AAB_H_

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <float.h>
#include <math.h>
#include <immintrin.h>
using namespace std;

namespace hispeed{

typedef struct range{
public:
	float mindist = 0;
	float maxdist = 0;

	bool operator>(range &d){
		return mindist>d.maxdist;
	}
	bool operator<(range &d){
		return maxdist<d.mindist;
	}
	bool operator==(range &d){
		return !(mindist>d.maxdist||maxdist<d.mindist);
	}
	friend std::ostream&
		operator<<(std::ostream& os, const range &d){
		os<<d.mindist<<"->"<<d.maxdist;
		return os;
	}
	void update(range &r){
		mindist = min(mindist, r.mindist);
		maxdist = min(maxdist, r.maxdist);
	}
	void print(){
		printf("[%f,%f]\n",mindist,maxdist);
	}
}range;

class aab{
public:
	float min[3];
	float max[3];

public:
	aab(){
		for(int i=0;i<3;i++){
			min[i] = DBL_MAX;
			max[i] = -DBL_MAX;
		}
	}

	aab(const aab &b){
		for(int i=0;i<3;i++){
			min[i] = b.min[i];
			max[i] = b.max[i];
		}
	}
	aab(float min_x, float min_y, float min_z,
			float max_x, float max_y, float max_z){
		min[0] = min_x;
		min[1] = min_y;
		min[2] = min_z;
		max[0] = max_x;
		max[1] = max_y;
		max[2] = max_z;
	}
	void set_box(float l0, float l1, float l2, float h0, float h1, float h2){
		min[0] = l0;
		min[1] = l1;
		min[2] = l2;
		max[0] = h0;
		max[1] = h1;
		max[2] = h2;
	}

	void update(float x, float y, float z){
		if(min[0]>x){
			min[0] = x;
		}
		if(max[0]<x){
			max[0] = x;
		}
		if(min[1]>y){
			min[1] = y;
		}
		if(max[1]<y){
			max[1] = y;
		}
		if(min[2]>z){
			min[2] = z;
		}
		if(max[2]<z){
			max[2] = z;
		}
	}

	void update(const aab &p){
		for(int i=0;i<3;i++){
			if(min[i]>p.min[i]){
				min[i] = p.min[i];
			}
			if(max[i]<p.max[i]){
				max[i] = p.max[i];
			}
		}
	}
	void set_box(const aab &b){
		for(int i=0;i<3;i++){
			min[i] = b.min[i];
			max[i] = b.max[i];
		}
	}


	bool intersect(aab &object){
		return !(object.min[0] >= max[0] || object.max[0] <= min[0] ||
				 object.min[1] >= max[1] || object.max[1] <= min[1] ||
				 object.min[2] >= max[2] || object.max[2] <= min[2]);
	}
	bool contains(aab *object){
		for(int i=0;i<3;i++){
			if(object->min[i]<min[i]){
				return false;
			}
			if(object->max[i]>max[i]){
				return false;
			}
		}
		return true;
	}

	bool contains(float *point){
		for(int i=0;i<3;i++){
			if(point[i]<min[i]){
				return false;
			}
			if(point[i]>max[i]){
				return false;
			}
		}
		return true;
	}

	friend std::ostream&
	operator<<(std::ostream& os, const aab &p){
		for(int i=0;i<3;i++){
		}
		os<<"(";
		os<<p.min[0]<<",";
		os<<p.min[1]<<",";
		os<<p.min[2]<<")";
		os<<" -> (";
		os<<p.max[0]<<",";
		os<<p.max[1]<<",";
		os<<p.max[2]<<")";
		return os;
	}
	float diagonal_length(){
		float dl = 0;
		for(int i=0;i<3;i++){
			dl += (max[i]-min[i])*(max[i]-min[i]);
		}
		return dl;
	}
	float volume(){
		return (max[0]-min[0])*(max[1]-min[1])*(max[2]-min[2]);
	}

	// get the possible minimum and maximum distance of
	// objects with their aabs
	range distance(const aab &b){
		range ret;
		ret.maxdist = 0;
		ret.mindist = 0;
		float tmp1 = 0;
		float tmp2 = 0;
		for(int i=0;i<3;i++){
			tmp1 = min[i]-b.max[i];
			tmp2 = max[i]-b.min[i];
			ret.maxdist += (tmp1+tmp2)*(tmp1+tmp2)/4;
			if(tmp2<0){
				ret.mindist += tmp2*tmp2;
			}else if(tmp1>0){
				ret.mindist += tmp1*tmp1;
			}
		}
		return ret;
	}
};


class weighted_aab:public aab{
public:
	int id;
	uint size = 1;
};


/*
 * each voxel contains the minimum boundary box
 * of a set of edges or triangles. It is an extension of
 * AAB with additional elements
 * */
class Voxel: public aab{
public:
	// point which the segments close with
	float core[3];
	// the pointer and size of the segment/triangle data in this voxel
	map<int, float *> data;
	map<int, int> size;
public:
	~Voxel(){
		reset();
	}
	void reset(){
		for(map<int, float *>::iterator it=data.begin();it!=data.end();it++){
			if(it->second!=NULL){
				delete []it->second;
				it->second = NULL;
			}
		}
		data.clear();
		size.clear();
	}
};


}


#endif /* HISPEED_AAB_H_ */
