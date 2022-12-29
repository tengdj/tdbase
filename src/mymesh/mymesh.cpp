/*
 * mymesh.cpp
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */

#include "mymesh.h"
#include <sstream>
#include <queue>

namespace hispeed{


Polyhedron::~Polyhedron(){

	for(Face *f:faces){
		delete f;
	}
	for(Vertex *p:vertices){
		assert(p->half_edges.size()==0 && p->opposite_half_edges.size()==0);
		delete p;
	}
	vertices.clear();
	faces.clear();
}

Vertex *Polyhedron::get_vertex(int vseq){
	assert(vseq>=0&&vseq<vertices.size());
	for(Vertex *v:vertices){
		if(vseq--==0){
			return v;
		}
	}
	assert(false);
	return NULL;
}

Face *Polyhedron::add_face(vector<Vertex *> &vs){
	Face *f = new Face(vs);
	faces.insert(f);
	return f;
}

Face *Polyhedron::remove_vertex(Vertex *v){

	vector<Face *> face_removed;
	for(Half_Edge *h:v->half_edges){
		face_removed.push_back(h->face);
	}

	unordered_set<Half_Edge *> boundary;
	for(Face *f:face_removed){
		for(Half_Edge *h:f->half_edges){
			if(h->end_vertex != v && h->vertex != v && h->opposite){
				if(h->opposite){
					boundary.insert(h->opposite);
				}
			}
		}
		faces.erase(f);
		delete f;
	}
	face_removed.clear();

	vertices.erase(v);

	assert(v->half_edges.size()==0 && v->opposite_half_edges.size()==0);
	vector<Vertex *> ordered_boundary;

	Half_Edge *cur = (*boundary.begin());
	Vertex *head = cur->end_vertex;
	do{
		boundary.erase(cur);
		ordered_boundary.push_back(cur->end_vertex);
		for(auto h:boundary){
			if(h->end_vertex == cur->vertex){
				cur = h;
				break;
			}
		}
	}while(!boundary.empty());
	Face *nface = add_face(ordered_boundary);
	nface->added = true;
	ordered_boundary.clear();

	return nface;
}

void Polyhedron::reset_states(){
	for(Vertex *v:vertices){
		v->added = false;
		v->removable = true;
	}

	for(Face *f:faces){
		f->added = false;
	}
}

void Polyhedron::compress(){
	Vertex *seed = get_vertex(0);
	queue<Vertex *> wq;
	wq.push(seed);
	while(!wq.empty()){
		Vertex *v = wq.front();
		wq.pop();
		for(Half_Edge *he:v->half_edges){
			if(v->removable){
				he->end_vertex->removable = false;
			}
			if(!he->end_vertex->added){
				he->end_vertex->added = true;
				wq.push(he->end_vertex);
			}
		}
		if(v->is_removable()){
			vector<Half_Edge *> hes;
			for(Half_Edge *he:v->half_edges){
				hes.push_back(he);
			}
			assert(hes.size()>1);
			for(Half_Edge *he:hes){
				Face *f = he->face;
				assert(f);
				if(f->facet_degree()>3){
					Face *nf = f->split(v);
					if(nf){
						faces.insert(nf);
					}
				}
			}
			remove_vertex(v);
		}
	}
}

int Polyhedron::remove_orphan_vertices(){
	vector<Vertex *> orphan;
	for(auto v:vertices){
		if(v->half_edges.size()==0){
			assert(v->opposite_half_edges.size()==0);
			orphan.push_back(v);
		}
	}
	int ret = orphan.size();
	for(auto v:orphan){
		vertices.erase(v);
		delete v;
	}
	orphan.clear();
	return ret;
}

// create a new half edge, setup the opposite of this half edge if needed
Half_Edge::Half_Edge(Vertex *v1, Vertex *v2){
	vertex = v1;
	end_vertex = v2;
	vertex->half_edges.insert(this);
	end_vertex->opposite_half_edges.insert(this);

	// in case this is the second half edge
	for(Half_Edge *h:v2->half_edges){
		if(h->end_vertex == v1){
			if(h->opposite){
				printf("create half edge:\n");
				v1->print();
				v2->print();
				h->opposite->face->print_off();
			}
			assert(h->opposite==NULL);
			h->opposite = this;
			this->opposite = h;
		}
	}
}

Half_Edge::~Half_Edge(){
	// reset the opposite;
	if(opposite!=NULL){
		assert(opposite->opposite = this);
		opposite->opposite = NULL;
	}

	// detach from the vertices
	assert(vertex && end_vertex);
	assert(vertex->half_edges.find(this)!=vertex->half_edges.end());
	assert(end_vertex->opposite_half_edges.find(this)!=end_vertex->opposite_half_edges.end());
	vertex->half_edges.erase(this);
	end_vertex->opposite_half_edges.erase(this);
}

void Face::remove(Half_Edge *rh){
	half_edges.erase(rh);
	for(Half_Edge *h:half_edges){
		if(h->next==rh){
			h->next = NULL;
		}
	}
	delete rh;
}

Face *Face::split(Vertex *v){
	if(facet_degree()==3){
		return NULL;
	}

//	printf("splitting\n");
//	v->print();
//	this->print_off();
//	v->print();
//	print_off();
//	printf("\n");
	Half_Edge *first = NULL;
	for(Half_Edge *h:half_edges){
		if(h->vertex==v){
			first = h;
			break;
		}
	}
	assert(first);

	if(this->degree()==4 && first->next->end_vertex->degree()==2){
		return NULL;
	}

	Half_Edge *last = first;
	Half_Edge *second_last = NULL;
	while(last->end_vertex!=v){
		second_last = last;
		last = last->next;
	}
	assert(last->end_vertex == v);

	Vertex *last_vertex = last->vertex;
	Vertex *first_vertex = first->end_vertex;
	Half_Edge *first_next = first->next;
	remove(first);
	remove(last);
	Half_Edge *new_he = new Half_Edge(last_vertex, first_vertex);
	new_he->face = this;
	half_edges.insert(new_he);

	second_last->next = new_he;
	new_he->next = first_next;

	vector<Vertex *> tmpvertices;
	for(Vertex *tv:vertices){
		if(tv!=v){
			tmpvertices.push_back(tv);
		}
	}
	vertices.clear();
	vertices.insert(vertices.begin(), tmpvertices.begin(), tmpvertices.end());

	vector<Vertex *> newfc;
	newfc.push_back(v);
	newfc.push_back(first_vertex);
	newfc.push_back(last_vertex);
	Face *f = new Face(newfc);
	return f;
}

}

