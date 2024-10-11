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
#include <map>
#include <assert.h>
#include <cstdint>
using namespace std;

namespace tdbase{

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
	void print(){
		printf("[%f,%f]\n",mindist,maxdist);
	}
	bool valid(){
		return mindist<=maxdist;
	}
}range;

class aab{
public:
	float low[3];
	float high[3];
public:
	aab();
	aab(const aab &b);
	aab(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z);
	void set_box(float l0, float l1, float l2, float h0, float h1, float h2);
	void reset();
	void update(float x, float y, float z);
	void update(const aab &p);
	void set_box(const aab &b);
	bool intersect(aab &object);
	bool contains(aab *object);
	bool contains(float *point);
	void print();

	friend std::ostream&
	operator<<(std::ostream& os, const aab &p){
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
	float diagonal_length();
	float volume();
	range distance(const aab &b);

	inline float length(){return high[0] - low[0];};
	inline float width(){return high[1] - low[1];};
	inline float height(){return high[2] - low[2];};

};


class weighted_aab:public aab{
public:
	int id;
	uint32_t size = 1;
};

/*
 * each voxel contains the minimum boundary box
 * of a set of edges or triangles. It is an extension of
 * AAB with additional elements
 * */
class Voxel: public aab{
public:
	// point which the triangles close with
	float core[3];

	float *triangles = NULL;
	float *hausdorff = NULL;
	size_t num_triangles = 0;
	size_t capacity = 0;

	// offset and volume information
	map<int, size_t> offset_lod;
	map<int, size_t> volume_lod;

	bool owned = false;
public:
	~Voxel();
	void clear();
	void reserve(size_t size);

	void insert(const float *t, const float *h);
	void batch_load(const float *t, const float *h, size_t s);
	void external_load(float *t, float *h, size_t s);
	void print();

	float getHausdorffDistance(int offset);
	float getProxyHausdorffDistance(int offset);
};


}


#endif /* HISPEED_AAB_H_ */
