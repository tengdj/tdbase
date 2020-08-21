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

/*
 *
 * extract the skeleton of current polyhedron
 *
 * */
Skeleton *HiMesh::extract_skeleton(){
	//assert(i_decompPercentage==100&&this->b_jobCompleted &&
	//		"the skeleton can only be extracted on the fully decompressed polyhedron");
	if(skeletons.find(i_decompPercentage)!=skeletons.end()){
		return skeletons[i_decompPercentage];
	}
	struct timeval start = get_cur_time();
	Skeleton *skeleton  = new Skeleton();
	std::stringstream os;
	os << *this;
	logt("dump to stream", start);
	Triangle_mesh tmesh;
	os >> tmesh;
	logt("convert to triangle mesh", start);
	assert(CGAL::Polygon_mesh_processing::triangulate_faces(faces(tmesh), tmesh));
	for(boost::graph_traits<Triangle_mesh>::face_descriptor fit : faces(tmesh))
	    if (next(next(halfedge(fit, tmesh), tmesh), tmesh)
	        !=   prev(halfedge(fit, tmesh), tmesh))
	      std::cerr << "Error: non-triangular face left in mesh." << std::endl;

	if (!CGAL::is_triangle_mesh(tmesh)){
		log("Input geometry is not triangulated.");
		exit(-1);
	}
	logt("triangulate", start);

	try{
		Skeletonization mcs(tmesh);
		logt("initialize skeletonization", start);

		// 1. Contract the mesh by mean curvature flow.
		mcs.contract_geometry();
		logt("contract geometry", start);

		// 2. Collapse short edges and split bad triangles.
		mcs.collapse_edges();
		logt("collapse edges", start);

		mcs.split_faces();
		logt("split faces", start);

		// 3. Fix degenerate vertices.
		//mcs.detect_degeneracies();

		// Perform the above three steps in one iteration.
		mcs.contract();
		logt("contract", start);

		// Iteratively apply step 1 to 3 until convergence.
		mcs.contract_until_convergence();
		logt("contract until convergence", start);

		// Convert the contracted mesh into a curve skeleton and
		// get the correspondent surface points
		mcs.convert_to_skeleton(*skeleton);
		logt("convert to skeleton", start);

	}catch(std::exception &exc){
		log(exc.what());
		exit(-1);
	}
	skeletons[i_decompPercentage] = skeleton;
	return skeleton;
}

vector<Point> HiMesh::get_skeleton_points(int num_cores){
	vector<Point> ret;
	Skeleton *skeleton = extract_skeleton();
	if(!skeleton){
		return ret;
	}
	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
		auto p = (*skeleton)[v].point;
		ret.push_back(Point(p.x(),p.y(),p.z()));
	}

	int num_points = ret.size();
	int skeleton_sample_rate = (num_cores*100)/num_points;
	if(num_cores>num_points){
		num_cores = num_points;
		return ret;
	}
	bool *need_keep = new bool[num_points];
	for(int i=0;i<num_points;i++){
		need_keep[i] = false;
	}
	for(int i=0;i<num_cores;i++){
		int target = 0;
		int try_round = 0;
		bool assigned = true;
		while(try_round++<10){
			target = hispeed::get_rand_number(num_cores);
			if(!need_keep[target]){
				need_keep[target] = true;
				assigned = true;
				break;
			}
		}
		if(!assigned){
			for(int j=target;j<num_points;j++){
				if(!need_keep[j]){
					need_keep[j] = true;
					assigned = true;
					break;
				}
			}
		}
		if(!assigned){
			for(int j=target;j>=0;j--){
				if(!need_keep[j]){
					need_keep[j] = true;
					assigned = true;
					break;
				}
			}
		}
		assert(assigned);
	}
	int index = 0;
	for(auto it=ret.begin();it!=ret.end();){
		if(!need_keep[index++]){
			ret.erase(it);
		}else{
			it++;
		}
	}
	delete []need_keep;
	return ret;
}

/*
 * get the skeleton of the polyhedron, and then group the triangles
 * or edges around the points of the skeleton, and generate a list
 * of axis aligned boxes for those sets
 * */
vector<Voxel *> HiMesh::generate_voxels(int voxel_size){
	timeval start = hispeed::get_cur_time();
	vector<Voxel *> voxels;
	int lod = i_decompPercentage;
	if(size_of_vertices()<voxel_size*3){
		Voxel *v = new Voxel();
		v->box = aab(bbMin[0],bbMin[1],bbMin[2],bbMax[0],bbMax[1],bbMax[2]);
		v->size[lod] = size_of_edges();
		voxels.push_back(v);
		return voxels;
	}
	// sample the points of the skeleton with the calculated sample rate
	int num_cores = size_of_vertices()/voxel_size;
	assert(num_cores>3);
	// this step takes 99 percent of the computation load
	vector<Point> skeleton_points = get_skeleton_points(num_cores);
	for(int i=0;i<skeleton_points.size();i++){
		Voxel *v = new Voxel();
		v->core[0] = skeleton_points[i][0];
		v->core[1] = skeleton_points[i][1];
		v->core[2] = skeleton_points[i][2];
		v->size[lod] = 0;
		voxels.push_back(v);
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

	for(Vertex_const_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		Point p = vit->point();
		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<skeleton_points.size();j++){
			float cur_dist = distance(skeleton_points[j], p);
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		voxels[gid]->box.update(p.x(),p.y(),p.z());
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

	size_t num_of_data = 0;
	int  size_of_datum = 0;
	float *data_buffer = NULL;
	int lod = i_decompPercentage;
	// the voxel should not be filled
	if(voxels[0]->data.find(lod)!=voxels[0]->data.end()){
		return;
	}
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
}

HiMesh::HiMesh(char* data, long length):
		MyMesh(0, DECOMPRESSION_MODE_ID, 12, data, length, true){
}
HiMesh::HiMesh(char* data, long length, bool od):
		MyMesh(0, DECOMPRESSION_MODE_ID, 12, data, length, od){
}
void HiMesh::get_segments(){
	segments.clear();
	for(Edge_const_iterator eit = edges_begin(); eit!=edges_end(); ++eit){
		Point p1 = eit->vertex()->point();
		Point p2 = eit->opposite()->vertex()->point();
		if(p1!=p2){
			segments.push_back(Segment(p1, p2));
		}
	}
}

void HiMesh::to_wkt(){
	cout<<"POLYHEDRALSURFACE Z (";
	bool lfirst = true;
	for ( Facet_const_iterator fit = facets_begin(); fit != facets_end(); ++fit){
		if(lfirst){
			lfirst = false;
		}else{
			cout<<",";
		}
		cout<<"((";
		bool first = true;
		Halfedge_around_facet_const_circulator hit(fit->facet_begin()), end(hit);
		Point firstpoint;
		do {
			Point p = hit->vertex()->point();
			if(!first){
				cout<<",";
			}else{
				firstpoint = p;
			}
			first = false;
			cout<<p[0]<<" "<<p[1]<<" "<<p[2];
			// Write the current vertex id.
		} while(++hit != end);
		cout<<","<<firstpoint[0]<<" "<<firstpoint[1]<<" "<<firstpoint[2];
		cout<<"))";
	}
	cout<<")"<<endl;
}

SegTree *HiMesh::get_aabb_tree(){
	advance_to(100);
	if(aabb_tree==NULL){
		get_segments();
		aabb_tree = new SegTree(segments.begin(), segments.end());
		aabb_tree->accelerate_distance_queries();
	}
	return aabb_tree;
}

float HiMesh::get_volume() {
	Nef_polyhedron inputpoly;
	stringstream ss;
	ss<<*this;
	CGAL::OFF_to_nef_3(ss, inputpoly);
	// to check if the intersected object can be converted to polyhedron or not
	std::vector<Polyhedron> PList;

	// decompose non-convex volume to convex parts
	cout<<"teng1"<<endl;

	convex_decomposition_3(inputpoly);
	cout<<"teng2"<<endl;

	for(Volume_const_iterator ci = ++inputpoly.volumes_begin() ; ci != inputpoly.volumes_end(); ++ci) {
		if(ci->mark()) {
			cout<<"teng"<<endl;
			Polyhedron P;
			inputpoly.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
			PList.push_back(P);
		}
	}
	cout<<PList.size()<<endl;

	double total_volume = 0, hull_volume = 0;
	for(Polyhedron poly:PList)
	{
		std::vector<Point> L;
		for (Polyhedron::Vertex_const_iterator  it = poly.vertices_begin(); it != poly.vertices_end(); it++) {
			L.push_back(Point(it->point().x(), it->point().y(), it->point().z()));
		}
		Triangulation T(L.begin(), L.end());
		hull_volume = 0;
		for(Triangulation::Finite_cells_iterator it = T.finite_cells_begin(); it != T.finite_cells_end(); it++) {
			Tetrahedron tetr = T.tetrahedron(it);
			hull_volume += to_double(tetr.volume());
		}
		total_volume += hull_volume;
	}
	return total_volume;
}

TriangleTree *get_aabb_tree(Polyhedron *p){
	TriangleTree *tree = new TriangleTree(faces(*p).first, faces(*p).second, *p);
	tree->accelerate_distance_queries();
	return tree;
}

void HiMesh_Wrapper::fill_voxels(enum data_type seg_tri, bool release_mesh){
	pthread_mutex_lock(&lock);
	// mesh could be NULL when multiple threads compete
	if(mesh){
		mesh->fill_voxel(voxels, seg_tri);
		// filled the maximum LOD, release the mesh
		if(release_mesh){
			delete mesh;
			mesh = NULL;
		}
	}
	pthread_mutex_unlock(&lock);
}




}
