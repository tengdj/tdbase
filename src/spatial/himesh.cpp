/*
 * voxel_group.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */


#include "himesh.h"


namespace hispeed{

Polyhedron *HiMesh::to_polyhedron(){
	stringstream ss;
	ss<<*this;
	Polyhedron *poly = new Polyhedron();
	ss>>*poly;
	return poly;
}

Polyhedron *HiMesh::to_triangulated_polyhedron(){
	Polyhedron *poly = to_polyhedron();
	CGAL::Polygon_mesh_processing::triangulate_faces(*poly);
	return poly;
}

void HiMesh::advance_to(int lod){
	i_decompPercentage = lod;
	b_jobCompleted = false;
	completeOperation();
}


HiMesh::HiMesh(MyMesh *mesh):
		MyMesh(0, DECOMPRESSION_MODE_ID, mesh->p_data, mesh->dataOffset){
}
HiMesh::HiMesh(char *data, size_t dsize):
		MyMesh(0, DECOMPRESSION_MODE_ID, data, dsize){
}

HiMesh::~HiMesh(){
	//release_buffer();
	this->clear_aabb_tree();
}

aab HiMesh::get_box(){
	aab b;
	for(Vertex_const_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		Point p = vit->point();
		b.update(p.x(), p.y(), p.z());
	}
	return b;
}

string HiMesh::to_wkt(){
	std::stringstream ss;
	ss<<"POLYHEDRALSURFACE Z (";
	bool lfirst = true;
	for ( Facet_const_iterator fit = facets_begin(); fit != facets_end(); ++fit){
		if(lfirst){
			lfirst = false;
		}else{
			ss<<",";
		}
		ss<<"((";
		bool first = true;
		Halfedge_around_facet_const_circulator hit(fit->facet_begin()), end(hit);
		Point firstpoint;
		do {
			Point p = hit->vertex()->point();
			if(!first){
				ss<<",";
			}else{
				firstpoint = p;
			}
			first = false;
			ss<<p[0]<<" "<<p[1]<<" "<<p[2];
			// Write the current vertex id.
		} while(++hit != end);
		ss<<","<<firstpoint[0]<<" "<<firstpoint[1]<<" "<<firstpoint[2];
		ss<<"))";
	}
	ss<<")";
	return ss.str();
}

size_t HiMesh::size_of_edges(){
	return size_of_halfedges()/2;
}

float HiMesh::get_volume() {
	Polyhedron *poly = this->to_triangulated_polyhedron();
	float volume = CGAL::Polygon_mesh_processing::volume(*poly);
	delete poly;
	return volume;
}

bool HiMesh::has_same_vertices(){
	return true;
}

}
