/*
 * himesh_retrieve_elements.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */

#include "himesh.h"


namespace hispeed{


size_t HiMesh::fill_segments(float *segments){
	assert(segments);
	size_t size = size_of_edges();
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


size_t HiMesh::fill_triangles(float *triangles){
	size_t size = size_of_facets();
	assert(triangles);
	float *cur_S = triangles;
	int inserted = 0;
	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		Point p1 = f->halfedge()->vertex()->point();
		Point p2 = f->halfedge()->next()->vertex()->point();
		Point p3 = f->halfedge()->next()->next()->vertex()->point();
		for(int i=0;i<3;i++){
			*cur_S = p1[i];
			cur_S++;
		}
		for(int i=0;i<3;i++){
			*cur_S = p2[i];//-p1[i];
			cur_S++;
		}
		for(int i=0;i<3;i++){
			*cur_S = p3[i];//-p1[i];
			cur_S++;
		}
		inserted++;
	}
	assert(inserted==size);
	return size;
}

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


size_t HiMesh::fill_topology(unsigned short *&topology){
	size_t size = this->size_of_facets();
	if(!topology){
		topology = new unsigned short[size*3];
	}

	unsigned short *cur = topology;
	int inserted = 0;
	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		*cur = (unsigned short)f->halfedge()->vertex()->getId();
		cur++;
		*cur = (unsigned short)f->halfedge()->next()->vertex()->getId();
		cur++;
		*cur = (unsigned short)f->halfedge()->next()->next()->vertex()->getId();
		cur++;
	}
	assert(inserted==size);
	return size;

}

// assign each segment(0) or triangle(1) to the proper voxel
size_t HiMesh::fill_voxels(vector<Voxel *> &voxels, enum data_type seg_or_triangle){
	assert(voxels.size()>0);

	size_t num_of_data = 0;
	int  size_of_datum = 0;
	float *data_buffer = NULL;
	int lod = i_decompPercentage;
	// the voxel should not be filled
//	if(voxels[0]->data.find(lod)!=voxels[0]->data.end()){
//		return;
//	}
	if(seg_or_triangle==DT_Segment){
		num_of_data = size_of_edges();
		size_of_datum = 6;
		data_buffer = new float[num_of_data*size_of_datum];
		fill_segments(data_buffer);
	}else{
		num_of_data = size_of_facets();
		size_of_datum = 9;
		data_buffer = new float[num_of_data*size_of_datum];
		fill_triangles(data_buffer);
	}

	// for the special case only one voxel exist
	if(voxels.size()==1){
		voxels[0]->size[lod] = num_of_data;
		voxels[0]->data[lod] = new float[num_of_data*size_of_datum];
		memcpy(voxels[0]->data[lod],
			   data_buffer,
			   num_of_data*size_of_datum*sizeof(float));
		delete []data_buffer;
		return num_of_data;
	}

	// now reorganize the data with the voxel information given
	// assign each segment to a proper group
	// we tried voronoi graph, but for some reason it's
	// even slower than the brute force method
	int *groups = new int[num_of_data];
	int *group_count = new int[voxels.size()];
	for(int i=0;i<voxels.size();i++){
		voxels[i]->size[lod] = 0;
		voxels[i]->data[lod] = NULL;
		group_count[i] = 0;
	}
	for(int i=0;i<num_of_data;i++){
		// for both segment and triangle, we assign it with only the first
		// point
		float *p1 = data_buffer+i*size_of_datum;
		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<voxels.size();j++){
			float cur_dist = distance(voxels[j]->core, p1);
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		groups[i] = gid;
		group_count[gid]++;
	}


	for(int i=0;i<voxels.size();i++){
		if(group_count[i]>0){
			voxels[i]->data[lod] = new float[group_count[i]*size_of_datum];
		}
	}

	// copy the data to the proper position in the segment_buffer
	for(int i=0;i<num_of_data;i++){
		Voxel *v = voxels[groups[i]];
		memcpy((void *)(v->data[lod]+v->size[lod]*size_of_datum),
			   (void *)(data_buffer+i*size_of_datum),
			   size_of_datum*sizeof(float));
		v->size[lod]++;
	}

	delete groups;
	delete group_count;
	delete data_buffer;
	return num_of_data;
}


}
