/*
 * aab.h
 *
 *  Created on: Nov 14, 2019
 *      Author: teng
 *
 *  define of Axis-Aligned Box
 */

#ifndef HISPEED_AAB_H_
#define HISPEED_AAB_H_

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <float.h>
using namespace std;

namespace hispeed{

class aab{
public:
	float min[3];
	float max[3];
	uint weight = 0;
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
		min[2] = min_y;
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


	bool intersect(aab *object){
		return !(object->min[0] > max[0] || object->max[0] < min[0]
		|| object->min[1] > max[1] || object->max[1] < min[1]
		|| object->min[2] > max[2] || object->max[2] < min[2]);
	}

	friend std::ostream&
	operator<<(std::ostream& os, const aab &p){
		for(int i=0;i<3;i++){
			os<<p.min[i]<<" ";
		}
		os<<"-> ";
		for(int i=0;i<3;i++){
			os<<p.max[i]<<" ";
		}
		os<<" | "<<p.weight;
		return os;
	}


};


}


#endif /* HISPEED_AAB_H_ */
