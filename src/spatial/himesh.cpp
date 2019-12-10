/*
 * voxel_group.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */



#include "himesh.h"


namespace hispeed{

void HiMesh::release_buffer(){
	for(map<int, float *>::iterator it=segment_buffer.begin();
			it!=segment_buffer.end();it++){
		if(it->second!=NULL){
			delete it->second;
			it->second=NULL;
		}
	}
	segment_buffer.clear();
	for(map<int, float *>::iterator it=triangle_buffer.begin();
			it!=triangle_buffer.end();it++){
		if(it->second!=NULL){
			delete it->second;
			it->second=NULL;
		}
	}
	triangle_buffer.clear();
}

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
		log("Input geometry is not triangulated.");
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
		//mcs.detect_degeneracies();

		// Perform the above three steps in one iteration.
		mcs.contract();
		// Iteratively apply step 1 to 3 until convergence.
		mcs.contract_until_convergence();
		// Convert the contracted mesh into a curve skeleton and
		// get the correspondent surface points
		mcs.convert_to_skeleton(*skeleton);
	}catch(std::exception &exc){
		log(exc.what());
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
	int lod = i_decompPercentage;
	if(size_of_edges()<voxel_size*3){
		Voxel *v = new Voxel();
		v->box = aab(bbMin[0],bbMin[1],bbMin[2],bbMax[0],bbMax[1],bbMax[2]);
		v->size[lod] = size_of_edges();
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
			v->size[lod] = 0;
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
		voxels[0]->size[lod] = size_of_edges();
		return voxels;
	}

	for(Edge_const_iterator eit = edges_begin(); eit!=edges_end(); ++eit){
		Point p1 = eit->vertex()->point();
		Point p2 = eit->opposite()->vertex()->point();

		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<skeleton_points.size();j++){
			float cur_dist = distance(skeleton_points[j], p1);
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		voxels[gid]->box.update(p1[0],p1[1],p1[2]);
		voxels[gid]->box.update(p2[0],p2[1],p2[2]);
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
	assert(triangles);
	size_t size = size_of_facets();
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
			*cur_S = p2[i]-p1[i];
			cur_S++;
		}
		for(int i=0;i<3;i++){
			*cur_S = p3[i]-p1[i];
			cur_S++;
		}
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
}



// the function to generate the segments(0) or triangle(1) and
// assign each segment(0) or triangle(1) to the proper voxel
void HiMesh::fill_voxel(vector<Voxel *> &voxels, enum data_type seg_or_triangle){
	assert(voxels.size()>0);
	release_buffer();
	size_t num_of_data = 0;
	int  size_of_datum = 0;
	float *data_buffer = NULL;
	int lod = i_decompPercentage;
	// already filled
	if(voxels[0]->data.find(lod)!=voxels[0]->data.end()){
		return;
	}
	if(seg_or_triangle==DT_Segment){
		num_of_data = size_of_edges();
		size_of_datum = 6;
		if(segment_buffer.find(lod)!=segment_buffer.end()){
			data_buffer = segment_buffer[lod];
		}else{
			data_buffer = new float[num_of_data*size_of_datum];
			fill_segments(data_buffer);
			segment_buffer[lod] = data_buffer;
		}
	}else{
		num_of_data = size_of_facets();
		size_of_datum = 9;
		if(triangle_buffer.find(lod)!=triangle_buffer.end()){
			data_buffer = triangle_buffer[lod];
		}else{
			data_buffer = new float[num_of_data*size_of_datum];
			fill_triangles(data_buffer);
			segment_buffer[lod] = data_buffer;
		}
	}

	// for the special case only one voxel exist
	if(voxels.size()==1){
		voxels[0]->size[lod] = num_of_data;
		voxels[0]->data[lod] = new float[num_of_data*size_of_datum];
		memcpy(voxels[0]->data[lod],
			   data_buffer,
			   num_of_data*size_of_datum*sizeof(float));
		return;
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
		group_count[i];
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
}

HiMesh::HiMesh(const char* data, long length):
		MyMesh(0, DECOMPRESSION_MODE_ID, 12, true, data, length){
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

SegTree *HiMesh::get_aabb_tree(){
	list<Segment> segments = get_segments();
	SegTree *tree = new SegTree(segments.begin(), segments.end());
	//tree->accelerate_distance_queries();
	return tree;
}

TriangleTree *get_aabb_tree(Polyhedron *p){
	TriangleTree *tree = new TriangleTree(faces(*p).first, faces(*p).second, *p);
	tree->accelerate_distance_queries();
	return tree;
}

void HiMesh_Wrapper::fill_voxels(enum data_type seg_tri){
	assert(mesh);
	pthread_mutex_lock(&lock);
	mesh->fill_voxel(voxels, seg_tri);
	pthread_mutex_unlock(&lock);
}


}
