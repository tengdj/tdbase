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

	Polyhedron *poly = this->to_polyhedron();
	if (!CGAL::is_triangle_mesh(*poly)){
		CGAL::Polygon_mesh_processing::triangulate_faces(*poly);
		if(global_ctx.verbose>=1){
			logt("convert to triangle mesh", start);
		}
	}
	os << *poly;
	delete poly;
	if(global_ctx.verbose>=1){
		logt("dump to stream", start);
	}

	os >> tmesh;
	assert(CGAL::is_triangle_mesh(tmesh));

	if(global_ctx.verbose>=1){
		logt("triangulate", start);
	}
	CGAL::extract_mean_curvature_flow_skeleton(tmesh, *skeleton);

//	try{
//		Skeletonization mcs(tmesh);
//		if(global_ctx.verbose){
//			logt("initialize skeletonization", start);
//		}
//
//		// 1. Contract the mesh by mean curvature flow.
//		mcs.contract_geometry();
//		if(global_ctx.verbose){
//			logt("contract geometry", start);
//		}
//
//		// 2. Collapse short edges and split bad triangles.
//		mcs.collapse_edges();
//		if(global_ctx.verbose){
//			logt("collapse edges", start);
//		}
//
//		mcs.split_faces();
//		if(global_ctx.verbose){
//			logt("split faces", start);
//		}
//
//		// 3. Fix degenerate vertices.
//		//mcs.detect_degeneracies();
//
//		// Perform the above three steps in one iteration.
//		mcs.contract();
//		if(global_ctx.verbose){
//			logt("contract", start);
//		}
//
//		// Iteratively apply step 1 to 3 until convergence.
//		mcs.contract_until_convergence();
//		if(global_ctx.verbose){
//			logt("contract until convergence", start);
//		}
//
//		// Convert the contracted mesh into a curve skeleton and
//		// get the correspondent surface points
//		mcs.convert_to_skeleton(*skeleton);
//		if(global_ctx.verbose){
//			logt("convert to skeleton", start);
//		}
//
//	}catch(std::exception &exc){
//		log(exc.what());
//	}
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
	aab box = get_mbb();
	if(voxel_num<=1){
		Voxel *v = new Voxel();
		v->set_box(box);
		voxels.push_back(v);
		return voxels;
	}

	//log("%d %d ",size_of_vertices(), voxel_num);
	//assert(num_cores>3);
	// this step takes 99 percent of the computation load
	vector<Point> skeleton_points = get_skeleton_points(voxel_num);
	for(int i=0;i<skeleton_points.size();i++){
		Voxel *v = new Voxel();
		v->core[0] = skeleton_points[i][0];
		v->core[1] = skeleton_points[i][1];
		v->core[2] = skeleton_points[i][2];
		voxels.push_back(v);
	}

	// return one single box if less than 2 points are sampled
	if(voxels.size()==0){
		Voxel *v = new Voxel();
		v->set_box(box);
		voxels.push_back(v);
		return voxels;
	}else if(voxels.size()==1){
		voxels[0]->set_box(box);
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
			float cur_dist = 0;
			for(int i=0;i<3;i++){
				cur_dist += (skeleton_points[j][i]-p[i])*(skeleton_points[j][i]-p[i]);
			}
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		voxels[gid]->update(p1.x(),p1.y(),p1.z());
		voxels[gid]->update(p2.x(),p2.y(),p2.z());
		voxels[gid]->update(p3.x(),p3.y(),p3.z());

		// we do not need the facets in each voxel, such that we do not need to insert one facet
		voxels[gid]->num_triangles++;
	}

	// erase the one without any data in it
	int vs=voxels.size();
	for(int i=0;i<vs;){
		if(voxels[i]->num_triangles==0){
			voxels.erase(voxels.begin()+i);
			vs--;
		}else{
			// reset to 0 as we never truly inserted any facets
			voxels[i]->num_triangles = 0;
			i++;
		}
	}
	skeleton_points.clear();

	return voxels;
}

}
