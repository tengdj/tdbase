/*
 * CGAL_common.cpp
 *
 *  Created on: Nov 13, 2019
 *      Author: teng
 */

#include "spatial.h"
#include "../geometry/aab.h"

namespace hispeed{

/*
 * extended from https://doc.cgal.org/latest/Polyhedron/index.html
 * return a axis aligned box from the min and max points given
 *
 * */
Polyhedron *make_cube(aab box) {

    // appends a cube of size [0,1]^3 to the polyhedron P.
	Polyhedron *P = new Polyhedron();

	Polyhedron::Halfedge_handle h =
    		P->make_tetrahedron(Point(box.max[0], box.min[1], box.min[2]),
								Point(box.min[0], box.min[1], box.max[2]),
								Point(box.min[0], box.min[1], box.min[2]),
								Point(box.min[0], box.max[1], box.min[2]));
	Polyhedron::Halfedge_handle g = h->next()->opposite()->next();
    P->split_edge( h->next());
    P->split_edge( g->next());
    P->split_edge( g);
    h->next()->vertex()->point()     = Point(box.max[0], box.min[1], box.max[2]);
    g->next()->vertex()->point()     = Point(box.min[0], box.max[1], box.max[2]);
    g->opposite()->vertex()->point() = Point(box.max[0], box.max[1], box.min[2]);
    Polyhedron::Halfedge_handle f = P->split_facet(g->next(),
                                      g->next()->next()->next());
    Polyhedron::Halfedge_handle e = P->split_edge(f);
    e->vertex()->point() = Point(box.max[0], box.max[1], box.max[2]);
    P->split_facet( e, f->next()->next());
    CGAL_postcondition( P->is_valid());
    return P;
}

void write_polyhedron(Polyhedron *mesh, const char *path){
	ofstream myfile;
	myfile.open(path);
	myfile << *mesh;
	myfile.close();
}

void write_box(aab box, int id, string prefix){
	char path[256];
	sprintf(path, "%s.%d.off", prefix.c_str(), id);
	hispeed::write_polyhedron(hispeed::make_cube(box), path);
}

void write_box(aab box, const char *path){
	hispeed::write_polyhedron(hispeed::make_cube(box), path);
}

void write_polyhedron(Polyhedron *mesh, int id){
	char path[256];
	sprintf(path, "/gisdata/%d.off", id);
	write_polyhedron(mesh, path);
}

inline string read_off_stdin(){
	string input_line;
	getline(std::cin, input_line);
	string whole_mesh = input_line;
	if(input_line.find("|") == std::string::npos){
		whole_mesh += "\n";
		while(getline(std::cin, input_line)){
			whole_mesh += input_line+"\n";
			//printf("%s\n",input_line.c_str());
		}
//		whole_mesh += "|";
//		printf("%s\n",whole_mesh.c_str());
//		printf("%ld\n",whole_mesh.size());
//		printf("%s\n",input_line.c_str());
	}else{
		boost::replace_all(whole_mesh, "|", "\n");
	}
	return whole_mesh;
}


string polyhedron_to_wkt(Polyhedron *poly){
	stringstream ss;
	ss.precision(10);
	ss<<"POLYHEDRALSURFACE Z (";
	bool lfirst = true;
	for ( Polyhedron::Facet_iterator fit = poly->facets_begin(); fit != poly->facets_end(); ++fit){
		if(lfirst){
			lfirst = false;
		}else{
			ss<<",";
		}
		ss<<"((";
		bool first = true;
		Polyhedron::Halfedge_around_facet_const_circulator hit(fit->facet_begin()), end(hit);
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
	ss<<")"<<endl;
	return ss.str();
}

MyMesh *get_mesh(string input_line, bool complete_compress){
	if(input_line.size()==0){
		return NULL;
	}
	if(input_line.find("|") != std::string::npos){
		boost::replace_all(input_line, "|", "\n");
	}
	char *data = new char[input_line.size()];
	memcpy(data, input_line.c_str(), input_line.size());

	// Init the random number generator.
	MyMesh *mesh = new MyMesh(100,
				 COMPRESSION_MODE_ID, 12,
				 data, input_line.size(), false);

	if(complete_compress){
		mesh->completeOperation();
	}
	return mesh;
}


MyMesh *read_mesh(){
	string mesh_str = read_off_stdin();
	MyMesh *mesh = get_mesh(mesh_str);
	assert(mesh && "this function must return a valid mesh");
	return mesh;
}

MyMesh *read_off(char *path){
	string str = read_file(path);
	return get_mesh(str);
}

Polyhedron *read_off_polyhedron(char *path){
	string input = read_file(path);
	stringstream ss;
	ss<<input;
	Polyhedron *poly = new Polyhedron();
	ss>>*poly;
	return poly;
}


MyMesh *decompress_mesh(MyMesh *compressed, int lod, bool complete_operation){
	MyMesh *decompressed = new MyMesh(lod,
			 DECOMPRESSION_MODE_ID, 12,
			 compressed->p_data, compressed->dataOffset, true);
	if(complete_operation){
		decompressed->completeOperation();
	}
	return decompressed;
}

MyMesh *decompress_mesh(char *data, size_t length, bool complete_operation){
	MyMesh *mesh = new MyMesh(100,
			DECOMPRESSION_MODE_ID, 12,
			data, length, true);
	if(complete_operation){
		mesh->completeOperation();
	}
	return mesh;
}

Polyhedron *read_polyhedron(){
	string input = read_off_stdin();
	printf("%s\n",input.c_str());
	stringstream ss;
	ss<<input;
	Polyhedron *poly = new Polyhedron();
	ss>>*poly;
	return poly;
}

//float get_volume(Polyhedron *polyhedron) {
//	Nef_polyhedron inputpoly;
//	stringstream ss;
//	ss<<*polyhedron;
//	cout<<"teng0"<<endl;
//	CGAL::OFF_to_nef_3(ss, inputpoly);
//	// to check if the intersected object can be converted to polyhedron or not
//	std::vector<Polyhedron> PList;
//
//	// decompose non-convex volume to convex parts
//	cout<<"teng1"<<endl;
//
//	convex_decomposition_3(inputpoly);
//	cout<<"teng2"<<endl;
//
//	for(Volume_const_iterator ci = ++inputpoly.volumes_begin() ; ci != inputpoly.volumes_end(); ++ci) {
//		if(ci->mark()) {
//			cout<<"teng"<<endl;
//			Polyhedron P;
//			inputpoly.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
//			PList.push_back(P);
//		}
//	}
//	cout<<PList.size()<<endl;
//
//	double total_volume = 0, hull_volume = 0;
//	for(Polyhedron poly:PList)
//	{
//		std::vector<Point> L;
//		for (Polyhedron::Vertex_const_iterator  it = poly.vertices_begin(); it != poly.vertices_end(); it++) {
//			L.push_back(Point(it->point().x(), it->point().y(), it->point().z()));
//		}
//		Triangulation T(L.begin(), L.end());
//		hull_volume = 0;
//		for(Triangulation::Finite_cells_iterator it = T.finite_cells_begin(); it != T.finite_cells_end(); it++) {
//			Tetrahedron tetr = T.tetrahedron(it);
//			hull_volume += to_double(tetr.volume());
//		}
//		total_volume += hull_volume;
//	}
//	return total_volume;
//}

}



