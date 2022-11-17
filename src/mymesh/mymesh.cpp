/*
 * mymesh.cpp
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */

#include "mymesh.h"
#include <sstream>

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
	Face *f = new Face();
	Half_Edge *prev = NULL;
	Half_Edge *head = NULL;
	for(int i=0;i<vs.size();i++){
		f->vertices.push_back(vs[i]);
		Vertex *nextv = vs[(i+1)%vs.size()];
		Half_Edge *hf = new Half_Edge(vs[i], nextv);
		hf->face = f;
		if(prev != NULL){
			prev->next = hf;
		}else{
			head = hf;
		}
		if(i==vs.size()-1){
			hf->next = head;
		}
		prev = hf;
		f->half_edges.insert(hf);
	}

	faces.insert(f);
	return f;
}

void Polyhedron::remove_vertex(Vertex *v){

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
	add_face(ordered_boundary);
	ordered_boundary.clear();
	return;
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

Half_Edge::Half_Edge(Vertex *v1, Vertex *v2){
	vertex = v1;
	end_vertex = v2;
	vertex->half_edges.insert(this);
	end_vertex->opposite_half_edges.insert(this);

	// in case this is the second half edge
	for(Half_Edge *h:v2->half_edges){
		if(h->end_vertex == v1){
			assert(h->opposite==NULL);
			h->opposite = this;
			this->opposite = h;
		}
	}
}

Half_Edge::~Half_Edge(){
	if(opposite!=NULL){
		assert(opposite->opposite = this);
		opposite->opposite = NULL;
	}
	assert(vertex && end_vertex);
	assert(vertex->half_edges.find(this)!=vertex->half_edges.end());
	vertex->half_edges.erase(this);
	assert(end_vertex->opposite_half_edges.find(this)!=end_vertex->opposite_half_edges.end());
	end_vertex->opposite_half_edges.erase(this);
}

void Face::remove(Half_Edge *rh){
	half_edges.erase(rh);
	for(Half_Edge *h:half_edges){
		if(h->next==rh){
			h->next = NULL;
		}
	}
}

Face *Face::split(Vertex *v){
	if(degree()==3){
		return NULL;
	}
	Half_Edge *first = NULL;
	for(Half_Edge *h:half_edges){
		if(h->vertex==v){
			first = h;
			break;
		}
	}
	assert(first);

	Half_Edge *last = first;
	while(last->end_vertex!=v){
		last = first->next;
	}
	assert(last->end_vertex == v);


}

}

