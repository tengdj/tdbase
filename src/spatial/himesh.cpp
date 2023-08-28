/*
 * voxel_group.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */


#include "himesh.h"


namespace hispeed{

int replacing_group::counter = 0;
int replacing_group::alive = 0;

HiMesh::HiMesh(string &str, bool completeop):
		CGAL::Polyhedron_3< CGAL::Simple_cartesian<float>, MyItems >(){

	struct timeval start = get_cur_time();
	boost::replace_all(str, "|", "\n");
	assert(str.size()!=0 && "input string should not be empty!");
	srand(PPMC_RANDOM_CONSTANT);
	i_mode = COMPRESSION_MODE_ID;
	// Create the compressed data buffer.
	const size_t d_capacity = 3*str.size();
	p_data = new char[d_capacity];
	// Fill the buffer with 0.
	std::istringstream is;
	is.str(str.c_str());
	is >> *this;
	if(size_of_facets()==0){
		std::cerr<<"failed to parse the OFF file into Polyhedron"<<endl;
		exit(EXIT_FAILURE);
	}

	if (keep_largest_connected_components(1) != 0){
		std::cerr << "Can't compress the mesh." << std::endl;
		std::cerr << "The codec doesn't handle meshes with several connected components." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!is_closed()){
		std::cerr << "Can't compress the mesh." << std::endl;
		std::cerr << "The codec doesn't handle meshes with borders." << std::endl;
		exit(EXIT_FAILURE);
	}

	compute_mbb();
	// Set the vertices of the edge that is the departure of the coding and decoding conquests.
	vh_departureConquest[0] = halfedges_begin()->opposite()->vertex();
	vh_departureConquest[1] = halfedges_begin()->vertex();
	//logt("load mesh", start);

	if(HiMesh::calculate_method == HCT_BVHTREE){
		get_aabb_tree_triangle();
		//logt("aabb", start);
	}

	if(HiMesh::calculate_method != HCT_NULL){
		updateVFMap();
		//logt("init triangles", start);
	}

	if(completeop){
		encode(0);
		//logt("encode", start);
	}
}

void HiMesh::updateVFMap(){

	//cout<<VFmap.size()<<endl;
	for(auto &vf:VFmap){
		vf.second.clear();
	}
	VFmap.clear();
	for(MyTriangle *tri:original_facets){
		delete tri;
	}
	original_facets.clear();

	for(Facet_iterator fit=facets_begin();fit!=facets_end();fit++){
		fit->tri = NULL;
	}

	const float area_unit = get_sample_density();
	for(Vertex_iterator vit=vertices_begin();vit!=vertices_end();vit++){
		Point p = vit->point();
		Halfedge_handle startH = vit->halfedge();
		Halfedge_handle h = startH->opposite(), end(h);

		vector<MyTriangle *> triangles;
		do {
			Face_handle f = h->face();
			if(f->tri == NULL){
				Triangle t(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
				f->tri = new MyTriangle(t, original_facets.size());
				original_facets.push_back(f->tri);
				sample_points(t, f->tri->sampled_points, area_unit);
			}
			triangles.push_back(f->tri);
		} while((h=h->opposite()->next()) != end);
		if(triangles.size()>0){
			VFmap[p] = triangles;
		}
	}
}

// in decompression mode
HiMesh::HiMesh(char *data, size_t dsize):
		CGAL::Polyhedron_3< CGAL::Simple_cartesian<float>, MyItems >(){
	assert(dsize>0);
	srand(PPMC_RANDOM_CONSTANT);

	i_mode = DECOMPRESSION_MODE_ID;
	p_data = new char[dsize];
	memcpy(p_data, data, dsize);
    readBaseMesh();
	// Set the vertices of the edge that is the departure of the coding and decoding conquests.
	vh_departureConquest[0] = vertices_begin();
	vh_departureConquest[1] = ++vertices_begin();
}

HiMesh::~HiMesh(){
	if(p_data!=NULL){
	   delete[] p_data;
	}
	clear_aabb_tree();
	for(replacing_group *rg:map_group){
		delete rg;
	}
	map_group.clear();
}

void HiMesh::compute_mbb(){
	mbb.reset();
	for(Vertex_const_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		Point p = vit->point();
		mbb.update(p.x(), p.y(), p.z());
	}
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

size_t HiMesh::size_of_triangles(){
	size_t tri_num = 0;
	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		tri_num += f->facet_degree()-2;
	}
	return tri_num;
}

Polyhedron *HiMesh::to_triangulated_polyhedron(){
	Polyhedron *poly = to_polyhedron();
	CGAL::Polygon_mesh_processing::triangulate_faces(*poly);
	return poly;
}

Polyhedron *HiMesh::to_polyhedron(){
	stringstream ss;
	ss<<*this;
	Polyhedron *poly = new Polyhedron();
	ss>>*poly;
	return poly;
}

HiMesh *HiMesh::clone_mesh(){
	stringstream ss;
	ss<<*this;
	string str = ss.str();
	HiMesh *nmesh = new HiMesh(str);
	return nmesh;
}

string HiMesh::to_off(){
	std::stringstream os;
	os << *this;
	return os.str();
}

void HiMesh::write_to_off(const char *path){
	string ct = to_off();
	hispeed::write_file(ct, path);
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

void HiMesh::write_to_wkt(const char *path){
	string ct = to_wkt();
	hispeed::write_file(ct, path);
}

/*
 * some other himesh common functions
 * */
aab HiMesh::shift(float x, float y, float z){
	for(Vertex_iterator vi=vertices_begin();vi!=vertices_end();vi++){
		Point p = vi->point();
		vi->point() = Point(p[0]+x, p[1]+y, p[2]+z);
	}
	compute_mbb();
	updateVFMap();
	return mbb;
}

aab HiMesh::shrink(float shrink){
	float mean_loc[3] = {(mbb.high[0]+mbb.low[0])/2,(mbb.high[1]+mbb.low[1])/2,(mbb.high[2]+mbb.low[2])/2};
	for(Vertex_iterator vi=vertices_begin();vi!=vertices_end();vi++){
		Point p = vi->point();
		vi->point() = Point((p[0]-mean_loc[0])/shrink+mean_loc[0],
							(p[1]-mean_loc[1])/shrink+mean_loc[1],
							(p[2]-mean_loc[2])/shrink+mean_loc[2]);

	}
	compute_mbb();
	updateVFMap();
	return mbb;
}

float HiMesh::area(){
	list<Triangle> triangles = this->get_triangles();
	float a = 0.0;
	for(const Triangle &t:triangles){
		a += hispeed::triangle_area(t);
	}
	return a;
}

}
