/*
 * CGAL_common.cpp
 *
 *  Created on: Nov 13, 2019
 *      Author: teng
 */

#include "spatial.h"
#include "../geometry/aab.h"

namespace hispeed{

/*
 * extended from https://doc.cgal.org/latest/Polyhedron/index.html
 * return a axis aligned box from the min and max points given
 *
 * */
Polyhedron *make_cube(aab box) {

    // appends a cube of size [0,1]^3 to the polyhedron P.
	Polyhedron *P = new Polyhedron();

	Polyhedron::Halfedge_handle h =
    		P->make_tetrahedron(Point(box.max[0], box.min[1], box.min[2]),
								Point(box.min[0], box.min[1], box.max[2]),
								Point(box.min[0], box.min[1], box.min[2]),
								Point(box.min[0], box.max[1], box.min[2]));
	Polyhedron::Halfedge_handle g = h->next()->opposite()->next();
    P->split_edge( h->next());
    P->split_edge( g->next());
    P->split_edge( g);
    h->next()->vertex()->point()     = Point(box.max[0], box.min[1], box.max[2]);
    g->next()->vertex()->point()     = Point(box.min[0], box.max[1], box.max[2]);
    g->opposite()->vertex()->point() = Point(box.max[0], box.max[1], box.min[2]);
    Polyhedron::Halfedge_handle f = P->split_facet(g->next(),
                                      g->next()->next()->next());
    Polyhedron::Halfedge_handle e = P->split_edge(f);
    e->vertex()->point() = Point(box.max[0], box.max[1], box.max[2]);
    P->split_facet( e, f->next()->next());
    CGAL_postcondition( P->is_valid());
    return P;
}

void write_polyhedron(Polyhedron *mesh, const char *path){
	ofstream myfile;
	myfile.open(path);
	myfile << *mesh;
	myfile.close();
}

void write_polyhedron(Polyhedron *mesh, int id){
	char path[256];
	sprintf(path, "%d.off", id);
	write_polyhedron(mesh, path);
}

MyMesh *get_mesh(string input_line, bool complete_compress){
	if(input_line.size()==0){
		return NULL;
	}
	boost::replace_all(input_line, "|", "\n");
	// Init the random number generator.
	MyMesh *mesh = new MyMesh(100,
				 COMPRESSION_MODE_ID, 12, true,
				 input_line.c_str(), input_line.size());
	if(complete_compress){
		mesh->completeOperation();
	}
	return mesh;
}


MyMesh *read_mesh(){
	string input_line;
	getline(std::cin, input_line);
	return get_mesh(input_line);
}


MyMesh *decompress_mesh(MyMesh *compressed, int lod){
	MyMesh *decompressed = new MyMesh(lod,
			 DECOMPRESSION_MODE_ID, 12, true,
			 compressed->p_data, compressed->dataOffset);
	decompressed->completeOperation();
	return decompressed;
}

MyMesh *decompress_mesh(char *data, size_t length, bool complete_operation){
	MyMesh *mesh = new MyMesh(100,
			DECOMPRESSION_MODE_ID, 12, true,
			data, length);
	if(complete_operation){
		mesh->completeOperation();
	}
	return mesh;
}

Polyhedron *read_polyhedron(){
	string input = read_polyhedron_str();
	stringstream ss;
	ss<<input;
	Polyhedron *poly = new Polyhedron();
	ss>>*poly;
	return poly;
}

}



