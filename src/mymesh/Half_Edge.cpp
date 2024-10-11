#include "mymesh.h"

namespace tmesh{
// create a new half edge, setup the opposite of this half edge if needed
Half_Edge::Half_Edge(Vertex *v1, Vertex *v2){
	vertex = v1;
	end_vertex = v2;
	vertex->half_edges.insert(this);
	end_vertex->opposite_half_edges.insert(this);

	// in case this is the second half edge
	for(Half_Edge *h:v2->half_edges){
		if(h->end_vertex == v1){
//			if(h->opposite){
//				printf("create half edge:\n");
//				v1->print();
//				v2->print();
//				h->opposite->face->print_off();
//			}
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

}
