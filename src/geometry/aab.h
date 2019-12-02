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
using namespace std;

namespace hispeed{

typedef struct range{
public:
	float closest = 0;
	float farthest = 0;
	bool operator>(range &d){
		return closest>d.farthest;
	}
	bool operator<(range &d){
		return farthest<d.closest;
	}
	bool operator==(range &d){
		return !(closest>d.farthest||farthest<d.closest);
	}
	friend std::ostream&
		operator<<(std::ostream& os, const range &d){
		os<<d.closest<<"->"<<d.farthest;
		return os;
	}
	void update(range &r){
		closest = max(closest, r.closest);
		farthest = min(farthest, r.farthest);
	}
}range;

class aab{
public:
	float min[3];
	float max[3];
	aab(){
		for(int i=0;i<3;i++){
			min[i] = DBL_MAX;
			max[i] = -DBL_MAX;
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
			max[1] = x;
		}
		if(min[2]>z){
			min[2] = z;
		}
		if(max[2]<z){
			max[2] = z;
		}
	}

	void update(aab &p){
		for(int i=0;i<3;i++){
			if(min[i]>p.min[i]){
				min[i] = p.min[i];
			}
			if(max[i]<p.max[i]){
				max[i] = p.max[i];
			}
		}
	}


	bool intersect(aab &object){
		return !(object.min[0] > max[0] || object.max[0] < min[0] ||
				 object.min[1] > max[1] || object.max[1] < min[1] ||
				 object.min[2] > max[2] || object.max[2] < min[2]);
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

	// get the possible minimum and maximum distance of
	// objects with their aabs
	range distance(const aab b){
		range ret;
		for(int i=0;i<3;i++){
			if(min[i]>b.max[i]){
				ret.closest += (min[i]-b.max[i])*(min[i]-b.max[i]);
			}else if(b.min[i]>max[i]){
				ret.closest += (b.min[i]-max[i])*(b.min[i]-max[i]);
			}
			float tmp = (b.max[i]+b.min[i]-max[i]-min[i])/2;
			ret.farthest += tmp*tmp;
			// else intersect in this dimension
		}
		return ret;
	}
};


class weighted_aab{
public:
	aab box;
	uint size;
};

}


#endif /* HISPEED_AAB_H_ */
