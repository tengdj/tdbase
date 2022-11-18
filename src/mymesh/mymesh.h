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

class Vertex;
class Face;
class Half_Edge;

class Vertex: public Point{
public:
	Vertex(float v1, float v2, float v3):Point(v1, v2, v3){}
	unordered_set<Half_Edge *> half_edges;
	unordered_set<Half_Edge *> opposite_half_edges;
	int id = 0;
	void print(){
		printf("%f %f %f\n", v[0], v[1], v[2]);
	}
};

class Half_Edge{
public:
	Vertex *vertex = NULL;
	Vertex *end_vertex = NULL;
	Face *face = NULL;
	Half_Edge *next = NULL;
	Half_Edge *opposite = NULL;
	Half_Edge(Vertex *v1, Vertex *v2);
	~Half_Edge();
};

class Face{
public:
	int id = 0;
	vector<Vertex *> vertices;
	unordered_set<Half_Edge *> half_edges;

	Face(){};
	~Face(){
		for(Half_Edge *h:half_edges){
			delete h;
		}
		half_edges.clear();
		vertices.clear();
	}

    Face(Vertex *v1, Vertex *v2, Vertex *v3){
    	vertices.push_back(v1);
    	vertices.push_back(v2);
    	vertices.push_back(v3);
    }

    Face(vector<Vertex *> &vs){
    	Half_Edge *prev = NULL;
    	Half_Edge *head = NULL;
    	for(int i=0;i<vs.size();i++){
    		vertices.push_back(vs[i]);
    		Vertex *nextv = vs[(i+1)%vs.size()];
    		Half_Edge *hf = new Half_Edge(vs[i], nextv);
    		hf->face = this;
    		if(prev != NULL){
    			prev->next = hf;
    		}else{
    			head = hf;
    		}
    		if(i==vs.size()-1){
    			hf->next = head;
    		}
    		prev = hf;
    		half_edges.insert(hf);
    	}
    }

    void print(){
    	for(Vertex *v:vertices){
    		v->print();
    	}
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

    int degree(){
    	return vertices.size();
    }
    // split the face and make sure the one without v as the new
    Face *split(Vertex *v);
    void remove(Half_Edge *h);

};



class Polyhedron{

public:
	int id = 0;
	unordered_set<Vertex *> vertices;
	unordered_set<Face *>	faces;

public:
	Polyhedron(int i=0){id=i;}
	~Polyhedron();

	// I/O
	void load(string path);
	bool parse(string str);
	bool parse(const char*, size_t);
	void dumpto(string path);
	void print();
	string to_string();
	Vertex *get_vertex(int vseq=0);

	// element operating
	Face *add_face(vector<Vertex *> &vs);
	Face *remove_vertex(Vertex *v);

	// mesh fixing
	int remove_orphan_vertices();
	void merge_vertex();
	bool fill_holes();
	void remove_redundant();
	vector<Polyhedron *> depart();
	void evaluate();

};

}
#endif /* SRC_MYMESH_MYMESH_H_ */
