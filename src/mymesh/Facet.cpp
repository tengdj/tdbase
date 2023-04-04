/*
 * Facet.cpp
 *
 *  Created on: Mar 24, 2023
 *      Author: teng
 */


#include "mymesh.h"

namespace tmesh{

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

	if(this->facet_degree()==4 && first->next->end_vertex->degree()==2){
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
