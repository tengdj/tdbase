/*
 * mymesh.h
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */

#ifndef SRC_MYMESH_MYMESH_H_
#define SRC_MYMESH_MYMESH_H_

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <set>
#include <map>
#include <stack>
#include <stdlib.h>
#include <tuple>
#include "util.h"

using namespace std;

namespace hispeed{

class Point{
public:
    float v[3];
    Point(float x1, float x2, float x3){
    	v[0] = x1;
    	v[1] = x2;
    	v[2] = x3;
    }
    Point(Point *pt){
    	assert(pt);
    	for(int i=0;i<3;i++){
    		v[i] = pt->v[i];
    	}
    };
    Point(){};
};

class Face;

class Vertex: public Point{
public:
	vector<Face *> faces;
};


class Face{
public:
	vector<int> vertices;

	Face(vector<int> &vs){
		vertices.insert(vertices.begin(), vs.begin(), vs.end());
	}
	Face(){};

    Face(int v1, int v2, int v3){
    	vertices.push_back(v1);
    	vertices.push_back(v2);
    	vertices.push_back(v3);
    }

    void print(){
    	for(int v:vertices){
    		printf("%d ",v);
    	}
    	printf("\n");
    }


    bool equal(const Face& rhs) const {
    	if(vertices.size()!=rhs.vertices.size()){
    		return false;
    	}
    	for(int i=0;i<vertices.size();i++){
    		if(vertices[i]!=rhs.vertices[i]){
    			return false;
    		}
    	}
    	return true;
    }
    bool operator==(const Face& rhs) const
    {
        return this->equal(rhs);
    }


};



class Polyhedron{

public:
	int id = 0;
	vector<Point *> points;
	vector<Face *>	faces;
	map<pair<int,int>, vector<Face *>> edges;

public:
	Polyhedron(int i=0){id=i;}
	~Polyhedron();

	// I/O
	void load(string path);
	bool parse(const char*, size_t);
	void dumpto(string path);
	void print();

	void get_edge();
	void add_face(Face *f);

	// mesh fixing
	void merge_vertex();
	void fill_holes();
	void remove_redundant();
	vector<Polyhedron *> depart();
	void evaluate();

};

}
#endif /* SRC_MYMESH_MYMESH_H_ */
