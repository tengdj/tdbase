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
	assert(i_decompPercentage==100&&this->b_jobCompleted &&
			"the skeleton can only be extracted on the fully decompressed polyhedron");
	struct timeval start = get_cur_time();
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

vector<Point> HiMesh::get_skeleton_points(){
	vector<Point> ret;
	Skeleton *skeleton = extract_skeleton();
	if(!skeleton){
		return ret;
	}
	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
		auto p = (*skeleton)[v].point;
		ret.push_back(Point(p.x(),p.y(),p.z()));
	}
	delete skeleton;
	return ret;
}

/*
 * get the skeleton of the polyhedron, and then group the triangles
 * or edges around the points of the skeleton, and generate a list
 * of axis aligned boxes for those sets
 * */
vector<Voxel *> HiMesh::generate_voxels(int voxel_size){
	vector<Voxel *> voxels;
	if(size_of_edges()<voxel_size*3){
		Voxel *v = new Voxel();
		v->box = aab(bbMin[0],bbMin[1],bbMin[2],bbMax[0],bbMax[1],bbMax[2]);
		v->size = size_of_edges();
		voxels.push_back(v);
		return voxels;
	}
	// this step takes 99 percent of the computation load
	std::vector<Point> skeleton_points = get_skeleton_points();

	// sample the points of the skeleton with the calculated sample rate
	int num_points = skeleton_points.size();
	int skeleton_sample_rate =
			((size_of_edges()/voxel_size)*100)/num_points;
	for(int i=0;i<num_points;){
		if(!hispeed::get_rand_sample(skeleton_sample_rate)){
			skeleton_points.erase(skeleton_points.begin()+i);
			num_points--;
		}else{
			Voxel *v = new Voxel();
			v->core[0] = skeleton_points[i][0];
			v->core[1] = skeleton_points[i][1];
			v->core[2] = skeleton_points[i][2];
			voxels.push_back(v);
			i++;
		}
	}

	if(voxels.size()==0){
		Voxel *v = new Voxel();
		voxels.push_back(v);
	}
	// return one single box if less than 2 points are sampled
	if(voxels.size()==1){
		voxels[0]->box = aab(bbMin[0],bbMin[1],bbMin[2],bbMax[0],bbMax[1],bbMax[2]);
		voxels[0]->size = size_of_edges();
		return voxels;
	}

	//building the voronoi graph
	Delaunay voronoi(skeleton_points.begin(), skeleton_points.end());

	// here we use a map to retrieve index from points
	map<Point, int> point_map;
	for(int i=0;i<skeleton_points.size();i++){
		point_map[skeleton_points[i]] = i;
	}
	for(Edge_const_iterator eit = edges_begin(); eit!=edges_end(); ++eit){
		Point p1 = eit->vertex()->point();
		Point p2 = eit->opposite()->vertex()->point();
		Point np = voronoi.nearest_vertex(p1)->point();
		voxels[point_map[np]]->box.update(p1[0],p1[1],p1[2]);
		voxels[point_map[np]]->box.update(p2[0],p2[1],p2[2]);
		voxels[point_map[np]]->size++;
	}

	// erase the one without any data in it
	int vs=voxels.size();
	for(int i=0;i<vs;){
		if(voxels[i]->size==0){
			voxels.erase(voxels.begin()+i);
			vs--;
		}else{
			i++;
		}
	}
	skeleton_points.clear();
	return voxels;
}

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
}

// the function to generate the segments and
// assign each segment to the proper voxel
void HiMesh::fill_voxel(vector<Voxel *> &voxels){
	assert(voxels.size()>0);
	if(segment_buffer){
		delete segment_buffer;
	}
	segment_buffer = new float[size_of_edges()*6];
	if(voxels.size()==1){
		fill_segments(segment_buffer);
		voxels[0]->size = size_of_edges();
		voxels[0]->data = segment_buffer;
		return;
	}
	float *tmp_segment_buffer = new float[size_of_edges()*6];
	int *groups = new int[size_of_edges()];
	fill_segments(tmp_segment_buffer);
	for(Voxel *v:voxels){
		v->size = 0;
	}
	// now reorganize the data with the voxel information given
	// assign each segment to a proper group
	// we tried voronoi graph, but for some reason it's
	// even slower than the brute force method
	for(int i=0;i<size_of_edges();i++){
		float *p1 = tmp_segment_buffer+i*6;
		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<voxels.size();j++){
			float cur_dist = 0;
			for(int t=0;t<3;t++){
				cur_dist += (voxels[j]->core[t]-p1[t])*(voxels[j]->core[t]-p1[t]);
			}
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		groups[i] = gid;
		voxels[gid]->size++;
	}

	// copy the data to the proper position in the segment_buffer
	int offset = 0;
	for(Voxel *v:voxels){
		v->data = segment_buffer+offset;
		offset += 6*(v->size);
	}
	for(int i=0;i<size_of_edges();i++){
		float *p1 = tmp_segment_buffer+i*6;
		float *p2 = tmp_segment_buffer+i*6+3;
		memcpy((void *)voxels[groups[i]]->data, (void*)(tmp_segment_buffer+i*6), 6*sizeof(float));
		voxels[groups[i]]->data += 6;
	}

	// reset pointer
	for(Voxel *v:voxels){
		v->data -= 6*(v->size);
	}
	delete groups;
	delete tmp_segment_buffer;
}

HiMesh::HiMesh(const char* data, long length):
		MyMesh(0, DECOMPRESSION_MODE_ID, 12, true, data, length){
}


}
