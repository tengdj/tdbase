/*
 * voxel_group.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */



#include "tile.h"


namespace hispeed{

/*
 *
 * extract the skeleton of current polyhedron
 *
 * */
Skeleton *HiMesh::extract_skeleton(){

	Skeleton *skeleton  = new Skeleton();
	std::stringstream os;
	os << *this;
	Triangle_mesh tmesh;
	os >> tmesh;
	if (!CGAL::is_triangle_mesh(tmesh)){
		std::cerr << "Input geometry is not triangulated." << std::endl;
		exit(-1);
	}
	try{
		Skeletonization mcs(tmesh);
		// 1. Contract the mesh by mean curvature flow.
		mcs.contract_geometry();
		// 2. Collapse short edges and split bad triangles.
		mcs.collapse_edges();
		mcs.split_faces();
		// 3. Fix degenerate vertices.
		mcs.detect_degeneracies();
		// Perform the above three steps in one iteration.
		mcs.contract();
		// Iteratively apply step 1 to 3 until convergence.
		mcs.contract_until_convergence();
		// Convert the contracted mesh into a curve skeleton and
		// get the correspondent surface points
		mcs.convert_to_skeleton(*skeleton);
	}catch(std::exception &exc){
		std::cerr<<exc.what()<<std::endl;
		exit(-1);
	}
	return skeleton;
}

/*
 * get the skeleton of the polyhedron, and then group the triangles
 * or edges around the points of the skeleton, and generate a list
 * of axis aligned boxes for those sets
 * */
std::vector<aab> HiMesh::generate_mbbs(int voxel_size){
	std::vector<aab> mbbs;
	if(size_of_facets()<voxel_size){
		mbbs.push_back(aab(bbMin[0],bbMin[1],bbMin[2],bbMax[0],bbMax[1],bbMax[2]));
		return mbbs;
	}
	std::vector<Point> skeleton_points;
	Skeleton *skeleton = extract_skeleton();
	assert(skeleton!=NULL);

	int num_points = 0;
	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
		num_points++;
	}
	// sample the points of the skeleton
	int skeleton_sample_rate =
			((size_of_facets()/voxel_size)*100)/num_points;


	// get the points in the skeleton
	// each mbb contains the same number of facets or edges
	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
		if(!hispeed::get_rand_sample(skeleton_sample_rate)){
			continue;
		}
		auto p = (*skeleton)[v].point;
		skeleton_points.push_back(Point(p.x(),p.y(),p.z()));
		mbbs.push_back(aab());
	}

	// reset the box_id of all facets
	for(Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
		fit->set_box_id(-1);
	}

	for(MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
		if(fit->get_box_id()==-1){
			Halfedge_handle heh = fit->halfedge();
			Halfedge_handle hIt = heh;
			do {
				Vertex_handle vh = hIt->vertex();
				Point p = vh->point();
				// update the corresponding mbb
				int min_index = 0;
				float min_dist = DBL_MAX;
				for(int i=0;i<skeleton_points.size();i++){
					float dist = get_distance(p,skeleton_points[i]);
					if(dist<min_dist){
						min_index = i;
						min_dist = dist;
					}
				}
				mbbs[min_index].update(p[0],p[1],p[2]);
				hIt = hIt->next();
			} while (hIt != heh);
		}
	}
	skeleton_points.clear();
	return mbbs;
}

size_t HiMesh::get_segments(float *segments){
	size_t size = size_of_halfedges()/2;
	if(segment_buffer==NULL){
		segment_buffer = new float[size*6*sizeof(float)];
		float *cur_S = segment_buffer;
		int inserted = 0;
		for(Edge_const_iterator eit = edges_begin(); eit!=edges_end(); ++eit){
			Point p1 = eit->vertex()->point();
			Point p2 = eit->opposite()->vertex()->point();
			assert(p1!=p2);
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
			*cur_S = p2.z();
			cur_S++;
			inserted++;
		}
		assert(inserted==size);
	}
	if(segments!=NULL){
		memcpy((void *)segments, (void *)segment_buffer, size*6*sizeof(float));
	}
	return size;
}

void HiMesh::advance_to(int lod){
	if(i_decompPercentage>=lod){
		return;
	}
	i_decompPercentage = lod;
	b_jobCompleted = false;
	completeOperation();
	// clean the buffer of segments if already assigned
	if(segment_buffer){
		delete segment_buffer;
		segment_buffer = NULL;
	}
	get_segments();
}

HiMesh::HiMesh(const char* data, long length):
		MyMesh(0, DECOMPRESSION_MODE_ID, 12, true, data, length){
}

}
