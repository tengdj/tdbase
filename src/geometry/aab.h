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
	bool operator>=(range &d){
		return mindist>=d.maxdist;
	}
	bool operator<(range &d){
		return maxdist<d.mindist;
	}
	bool operator<=(range &d){
		return maxdist<=d.mindist;
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
	float low[3];
	float high[3];
public:
	aab(){
		for(int i=0;i<3;i++){
			low[i] = DBL_MAX;
			high[i] = -DBL_MAX;
		}
	}

	aab(const aab &b){
		for(int i=0;i<3;i++){
			low[i] = b.low[i];
			high[i] = b.high[i];
		}
	}
	aab(float min_x, float min_y, float min_z,
			float max_x, float max_y, float max_z){
		low[0] = min_x;
		low[1] = min_y;
		low[2] = min_z;
		high[0] = max_x;
		high[1] = max_y;
		high[2] = max_z;
	}
	void set_box(float l0, float l1, float l2, float h0, float h1, float h2){
		low[0] = l0;
		low[1] = l1;
		low[2] = l2;
		high[0] = h0;
		high[1] = h1;
		high[2] = h2;
	}

	void update(float x, float y, float z){
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

	void update(const aab &p){
		for(int i=0;i<3;i++){
			if(low[i]>p.low[i]){
				low[i] = p.low[i];
			}
			if(high[i]<p.high[i]){
				high[i] = p.high[i];
			}
		}
	}
	void set_box(const aab &b){
		for(int i=0;i<3;i++){
			low[i] = b.low[i];
			high[i] = b.high[i];
		}
	}


	bool intersect(aab &object){
		return !(object.low[0] >= high[0] || object.high[0] <= low[0] ||
				 object.low[1] >= high[1] || object.high[1] <= low[1] ||
				 object.low[2] >= high[2] || object.high[2] <= low[2]);
	}
	bool contains(aab *object){
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

	bool contains(float *point){
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

	friend std::ostream&
	operator<<(std::ostream& os, const aab &p){
		for(int i=0;i<3;i++){
		}
		os<<"(";
		os<<p.low[0]<<",";
		os<<p.low[1]<<",";
		os<<p.low[2]<<")";
		os<<" -> (";
		os<<p.high[0]<<",";
		os<<p.high[1]<<",";
		os<<p.high[2]<<")";
		return os;
	}
	float diagonal_length(){
		float dl = 0;
		for(int i=0;i<3;i++){
			dl += (high[i]-low[i])*(high[i]-low[i]);
		}
		return dl;
	}
	float volume(){
		return (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
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
			tmp1 = low[i]-b.high[i];
			tmp2 = high[i]-b.low[i];
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
