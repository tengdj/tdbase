/*
 * himesh_internal_index.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */

#include "himesh.h"

namespace hispeed{

/*
 *
 * extract the skeleton of current polyhedron
 *
 * */
Skeleton *HiMesh::extract_skeleton(){

	struct timeval start = get_cur_time();
	Skeleton *skeleton  = new Skeleton();
	Triangle_mesh tmesh;
	std::stringstream os;

	//logt("dump to stream", start);
	//logt("convert to triangle mesh", start);
	if (!CGAL::is_triangle_mesh(tmesh)){
		os << *this;
	}else{
		Polyhedron *poly = to_triangulated_polyhedron();
		os << *poly;
		delete poly;
	}
	os >> tmesh;
	assert(CGAL::is_triangle_mesh(tmesh));

	//logt("triangulate", start);

	try{
		Skeletonization mcs(tmesh);
		//logt("initialize skeletonization", start);

		// 1. Contract the mesh by mean curvature flow.
		mcs.contract_geometry();
		//logt("contract geometry", start);

		// 2. Collapse short edges and split bad triangles.
		mcs.collapse_edges();
		//logt("collapse edges", start);

		mcs.split_faces();
		//logt("split faces", start);

		// 3. Fix degenerate vertices.
		//mcs.detect_degeneracies();

		// Perform the above three steps in one iteration.
		mcs.contract();
		//logt("contract", start);

		// Iteratively apply step 1 to 3 until convergence.
		mcs.contract_until_convergence();
		//logt("contract until convergence", start);

		// Convert the contracted mesh into a curve skeleton and
		// get the correspondent surface points
		mcs.convert_to_skeleton(*skeleton);
		//logt("convert to skeleton", start);

	}catch(std::exception &exc){
		log(exc.what());
		exit(-1);
	}
	return skeleton;
}


vector<Point> HiMesh::get_skeleton_points(int num_cores){
	vector<Point> ret;
	Skeleton *skeleton = extract_skeleton();
	if(!skeleton){
		return ret;
	}
	int num_vertices = 0;
	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
		num_vertices++;
	}
	double skeleton_sample_rate = (num_cores*1.0)/num_vertices;

	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
		if(tryluck(skeleton_sample_rate)){
			auto p = (*skeleton)[v].point;
			ret.push_back(Point(p.x(),p.y(),p.z()));
		}
	}
	return ret;
}

/*
 * get the skeleton of the polyhedron, and then group the triangles
 * or edges around the points of the skeleton, and generate a list
 * of axis aligned boxes for those sets
 * */
vector<Voxel *> HiMesh::generate_voxels_skeleton(int voxel_num){

	voxel_num = std::max(1, voxel_num);
	timeval start = hispeed::get_cur_time();
	vector<Voxel *> voxels;
	int lod = i_decompPercentage;
	aab box = get_box();
	if(voxel_num<=1){
		Voxel *v = new Voxel();
		v->set_box(box);
		v->size[lod] = size_of_edges();
		voxels.push_back(v);
		return voxels;
	}

	//log("%d %d %d cores",size_of_vertices(), voxel_size, num_cores);
	//assert(num_cores>3);
	// this step takes 99 percent of the computation load
	vector<Point> skeleton_points = get_skeleton_points(voxel_num);
	for(int i=0;i<skeleton_points.size();i++){
		Voxel *v = new Voxel();
		v->core[0] = skeleton_points[i][0];
		v->core[1] = skeleton_points[i][1];
		v->core[2] = skeleton_points[i][2];
		v->size[lod] = 0;
		voxels.push_back(v);
	}
	// return one single box if less than 2 points are sampled
	if(voxels.size()==0){
		Voxel *v = new Voxel();
		v->set_box(box);
		v->size[lod] = size_of_edges();
		voxels.push_back(v);
		return voxels;
	}else if(voxels.size()==1){
		voxels[0]->set_box(box);
		voxels[0]->size[lod] = size_of_edges();
		return voxels;
	}

	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		Point p1 = f->halfedge()->vertex()->point();
		Point p2 = f->halfedge()->next()->vertex()->point();
		Point p3 = f->halfedge()->next()->next()->vertex()->point();
		Point p((p1[0]+p2[0]+p3[0])/3, (p1[1]+p2[1]+p3[1])/3,(p1[2]+p2[2]+p3[2])/3);
		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<skeleton_points.size();j++){
			float cur_dist = distance(skeleton_points[j], p);
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		voxels[gid]->update(p1.x(),p1.y(),p1.z());
		voxels[gid]->update(p2.x(),p2.y(),p2.z());
		voxels[gid]->update(p3.x(),p3.y(),p3.z());
		voxels[gid]->size[lod]++;
	}

	// erase the one without any data in it
	int vs=voxels.size();
	for(int i=0;i<vs;){
		if(voxels[i]->size[lod]==0){
			voxels.erase(voxels.begin()+i);
			vs--;
		}else{
			i++;
		}
	}

	skeleton_points.clear();
	return voxels;
}

vector<Voxel *> HiMesh::voxelization(int voxel_size){
	vector<Voxel *> voxels;
	if(voxel_size<=1){
		Voxel *vox = new Voxel();
		vox->set_box(get_box());
		voxels.push_back(vox);
	}

	aab box = get_box();
	float min_dim = std::min(box.max[2]-box.min[2], std::min(box.max[1]-box.min[1], box.max[0]-box.min[0]));
	float div = (box.max[2]-box.min[2])*(box.max[1]-box.min[1])*(box.max[0]-box.min[0])/(min_dim*min_dim*min_dim);
	float multi = std::pow(1.0*voxel_size/div, 1.0/3);

	int dim[3];
	for(int i=0;i<3;i++){
		dim[i] = ((box.max[i]-box.min[i])*multi/min_dim+0.5);
		assert(dim[i]>0);
	}

	bool *taken = new bool[dim[0]*dim[1]*dim[2]];
	for(int i=0;i<dim[0]*dim[1]*dim[2];i++){
		taken[i] = false;
	}

//	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
//		Point p1 = f->halfedge()->vertex()->point();
//		Point p2 = f->halfedge()->next()->vertex()->point();
//		Point p3 = f->halfedge()->next()->next()->vertex()->point();
//	}
	for(Vertex_const_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		Point p = vit->point();
		int x = (p.x()-box.min[0])*dim[0]/(box.max[0]-box.min[0]);
		int y = (p.y()-box.min[1])*dim[1]/(box.max[1]-box.min[1]);
		int z = (p.z()-box.min[2])*dim[2]/(box.max[2]-box.min[2]);
		int idx = z*dim[1]*dim[0]+y*dim[0]+x;
		if(!taken[idx]){
			Voxel *vox = new Voxel();
			vox->min[0] = x*(box.max[0]-box.min[0])/dim[0];
			vox->min[1] = y*(box.max[1]-box.min[1])/dim[1];
			vox->min[2] = z*(box.max[2]-box.min[2])/dim[2];
			vox->max[0] = (x+1)*(box.max[0]-box.min[0])/dim[0];
			vox->max[1] = (y+1)*(box.max[1]-box.min[1])/dim[1];
			vox->max[2] = (z+1)*(box.max[2]-box.min[2])/dim[2];
			voxels.push_back(vox);
		}
		taken[idx] = true;
	}
	return voxels;
}




}
