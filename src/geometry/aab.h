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

typedef CGAL::Simple_cartesian<float> MyKernel;
typedef MyKernel::Point_3 Point;

namespace hispeed{

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
	aab(Point min, Point max){
		for(int i=0;i<3;i++){
			this->min[i] = min[i];
			this->max[i] = max[i];
		}
	}
	aab(std::vector<Point> &points){
		for(Point p:points){
			update(p);
		}
	}

	void update(Point &p){
		for(int i=0;i<3;i++){
			if(min[i]>p[i]){
				min[i] = p[i];
			}
			if(max[i]<p[i]){
				max[i] = p[i];
			}
		}
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
		return os;
	}

};


}


#endif /* HISPEED_AAB_H_ */
