/*
 * himesh_retrieve.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */

#include "himesh.h"


namespace tdbase{

/*
 * retrieve the surface elements
 *
 * */
size_t HiMesh::fill_vertices(float *&vertices){
	size_t size = this->size_of_vertices();
	if(!vertices){
		vertices = new float[size*3];
	}

	float *cur = vertices;

	for(Vertex_const_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		Point p = vit->point();
		*cur = p.x();
		cur++;
		*cur = p.y();
		cur++;
		*cur = p.z();
		cur++;
	}
	return size;
}

size_t HiMesh::fill_segments(float *&segments){
	size_t size = size_of_edges();
	segments = new float[size*6];
	float *cur_S = segments;
	int inserted = 0;
	for(Edge_const_iterator eit = edges_begin(); eit!=edges_end(); ++eit){
		Point p1 = eit->vertex()->point();
		Point p2 = eit->opposite()->vertex()->point();
		*cur_S = p1.x();
		cur_S++;
		*cur_S = p1.y();
		cur_S++;
		*cur_S = p1.z();
		cur_S++;
		*cur_S = p2.x();
		cur_S++;
		*cur_S = p2.y();
		cur_S++;
		*cur_S = p2.z()+0.1*(p1==p2); // pad one a little bit if the edge is a point
		cur_S++;
		inserted++;
	}
	assert(inserted==size);
	return size;
}

size_t HiMesh::fill_triangles(float *&triangles){
	size_t size = size_of_triangles();
	triangles = new float[9*size];
	assert(triangles);
	float *cur_S = triangles;
	int inserted = 0;
	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		Halfedge_const_handle e1 = f->halfedge();
		Halfedge_const_handle e2 = e1->next();
		do{
			Point p1 = e1->vertex()->point();
			Point p2 = e2->vertex()->point();
			Point p3 = e2->next()->vertex()->point();
			for(int i=0;i<3;i++){
				*cur_S = p1[i];
				cur_S++;
			}
			for(int i=0;i<3;i++){
				*cur_S = p2[i];
				cur_S++;
			}
			for(int i=0;i<3;i++){
				*cur_S = p3[i];
				cur_S++;
			}
			inserted++;
			e2 = e2->next();
		}while(e1!=e2->next());
	}
	assert(inserted==size);
	return size;
}

list<Point> HiMesh::get_vertices(){
	list<Point> vertices;
	for(HiMesh::Vertex_iterator v = vertices_begin(); v != vertices_end(); ++v){
		vertices.push_back(v->point());
	}
	return vertices;
}

list<Segment> HiMesh::get_segments(){
	list<Segment> segments;
	for(Edge_const_iterator eit = edges_begin(); eit!=edges_end(); ++eit){
		Point p1 = eit->vertex()->point();
		Point p2 = eit->opposite()->vertex()->point();
		if(p1!=p2){
			segments.push_back(Segment(p1, p2));
		}
	}
	return segments;
}

list<Triangle> HiMesh::get_triangles(){
	list<Triangle> triangles;
	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		Halfedge_const_handle e1 = f->halfedge();
		Halfedge_const_handle e2 = e1->next();
		do{
			triangles.push_back(Triangle(e1->vertex()->point(),
										 e2->vertex()->point(),
										 e2->next()->vertex()->point()));
			e2 = e2->next();
		}while(e1!=e2->next());
	}
	return triangles;
}

size_t HiMesh::fill_hausdorf_distances(float *&hausdorf){
	size_t size = size_of_triangles();
	hausdorf = new float[2*size];
	assert(hausdorf);
	float *cur_S = hausdorf;
	int inserted = 0;
	for ( Facet_iterator f = facets_begin(); f != facets_end(); ++f){
		Halfedge_const_handle e1 = f->halfedge();
		Halfedge_const_handle e2 = e1->next();
		do{
			*cur_S = f->getProxyHausdorff();
			cur_S++;
			*cur_S = f->getHausdorff();
			cur_S++;
			inserted++;
			e2 = e2->next();
		}while(e1!=e2->next());
	}
	assert(inserted==size);
	return size;
}

/*
 *
 * fill the triangles and their associated Hausdorff distances into proper voxels
 *
 * */
size_t HiMesh::fill_voxels(vector<Voxel *> &voxels){
	assert(voxels.size()>0);

	size_t num_of_element = size_of_triangles();
	float *triangle_buffer = NULL;
	float *hausdorf_buffer = NULL;
	fill_triangles(triangle_buffer);
	fill_hausdorf_distances(hausdorf_buffer);

	// for the special case only one voxel exist
	if(voxels.size()==1){
		voxels[0]->batch_load(triangle_buffer, hausdorf_buffer, num_of_element);
		delete []triangle_buffer;
		delete []hausdorf_buffer;
		return num_of_element;
	}

	// now reorganize the data with the voxel information given
	// assign each triangle to a proper group
	// we tried voronoi graph, but it's even slower than the brute force method
	int *groups = new int[num_of_element];
	int *group_count = new int[voxels.size()];
	for(int i=0;i<voxels.size();i++){
		group_count[i] = 0;
	}
	for(int i=0;i<num_of_element;i++){
		// for both segment and triangle, we assign it with only the first point
		float *p1 = triangle_buffer+i*9;
		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<voxels.size();j++){
			float cur_dist = 0;
			for(int t=0;t<3;t++){
				cur_dist += (p1[t]-voxels[j]->core[t])*(p1[t]-voxels[j]->core[t]);
			}
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		groups[i] = gid;
		group_count[gid]++;
	}

	// reserve enough space for holding the assigned triangles
	for(int i=0;i<voxels.size();i++){
		voxels[i]->reserve(group_count[i]);
	}

	// copy the data to the proper position in the buffers
	for(int i=0;i<num_of_element;i++){
		voxels[groups[i]]->insert(triangle_buffer+i*9, hausdorf_buffer+i*2);
	}

	// clean the buffer
	delete []groups;
	delete []group_count;
	delete []triangle_buffer;
	delete []hausdorf_buffer;
	return num_of_element;
}

}
