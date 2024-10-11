/*
 * CGAL_common.cpp
 *
 *  Created on: Nov 13, 2019
 *      Author: teng
 */

#include "himesh.h"
#include "aab.h"

namespace tdbase{

/*
 * extended from https://doc.cgal.org/latest/Polyhedron/index.html
 * return a axis aligned box from the min and max points given
 *
 * */
Polyhedron *make_cube(aab box) {

    // appends a cube of size [0,1]^3 to the polyhedron P.
	Polyhedron *P = new Polyhedron();

	Polyhedron::Halfedge_handle h =
    		P->make_tetrahedron(Point(box.high[0], box.low[1], box.low[2]),
								Point(box.low[0], box.low[1], box.high[2]),
								Point(box.low[0], box.low[1], box.low[2]),
								Point(box.low[0], box.high[1], box.low[2]));
	Polyhedron::Halfedge_handle g = h->next()->opposite()->next();
    P->split_edge( h->next());
    P->split_edge( g->next());
    P->split_edge( g);
    h->next()->vertex()->point()     = Point(box.high[0], box.low[1], box.high[2]);
    g->next()->vertex()->point()     = Point(box.low[0], box.high[1], box.high[2]);
    g->opposite()->vertex()->point() = Point(box.high[0], box.high[1], box.low[2]);
    Polyhedron::Halfedge_handle f = P->split_facet(g->next(),
                                      g->next()->next()->next());
    Polyhedron::Halfedge_handle e = P->split_edge(f);
    e->vertex()->point() = Point(box.high[0], box.high[1], box.high[2]);
    P->split_facet( e, f->next()->next());
    CGAL_postcondition( P->is_valid());
    return P;
}

Polyhedron *make_cubes(vector<aab *> &boxes) {

    // appends a cube of size [0,1]^3 to the polyhedron P.
	Polyhedron *P = new Polyhedron();

	for(aab *b:boxes){
		aab box = *b;
		Polyhedron::Halfedge_handle h =
	    		P->make_tetrahedron(Point(box.high[0], box.low[1], box.low[2]),
									Point(box.low[0], box.low[1], box.high[2]),
									Point(box.low[0], box.low[1], box.low[2]),
									Point(box.low[0], box.high[1], box.low[2]));
		Polyhedron::Halfedge_handle g = h->next()->opposite()->next();
	    P->split_edge( h->next());
	    P->split_edge( g->next());
	    P->split_edge( g);
	    h->next()->vertex()->point()     = Point(box.high[0], box.low[1], box.high[2]);
	    g->next()->vertex()->point()     = Point(box.low[0], box.high[1], box.high[2]);
	    g->opposite()->vertex()->point() = Point(box.high[0], box.high[1], box.low[2]);
	    Polyhedron::Halfedge_handle f = P->split_facet(g->next(),
	                                      g->next()->next()->next());
	    Polyhedron::Halfedge_handle e = P->split_edge(f);
	    e->vertex()->point() = Point(box.high[0], box.high[1], box.high[2]);
	    P->split_facet( e, f->next()->next());
	    CGAL_postcondition( P->is_valid());
	}
    return P;
}

Polyhedron *make_cubes(vector<Voxel *> &boxes) {
	 // appends a cube of size [0,1]^3 to the polyhedron P.
	Polyhedron *P = new Polyhedron();

	for(Voxel *b:boxes){
		Polyhedron::Halfedge_handle h =
				P->make_tetrahedron(Point(b->high[0], b->low[1], b->low[2]),
									Point(b->low[0], b->low[1], b->high[2]),
									Point(b->low[0], b->low[1], b->low[2]),
									Point(b->low[0], b->high[1], b->low[2]));
		Polyhedron::Halfedge_handle g = h->next()->opposite()->next();
		P->split_edge( h->next());
		P->split_edge( g->next());
		P->split_edge( g);
		h->next()->vertex()->point()     = Point(b->high[0], b->low[1], b->high[2]);
		g->next()->vertex()->point()     = Point(b->low[0], b->high[1], b->high[2]);
		g->opposite()->vertex()->point() = Point(b->high[0], b->high[1], b->low[2]);
		Polyhedron::Halfedge_handle f = P->split_facet(g->next(),
										  g->next()->next()->next());
		Polyhedron::Halfedge_handle e = P->split_edge(f);
		e->vertex()->point() = Point(b->high[0], b->high[1], b->high[2]);
		P->split_facet( e, f->next()->next());
		CGAL_postcondition( P->is_valid());
	}
	return P;
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

/*
 * persistent functions
 * */
void write_polyhedron(Polyhedron *mesh, const char *path){
	ofstream myfile;
	myfile.open(path);
	myfile << *mesh;
	myfile.close();
}

void write_box(aab box, int id, string prefix){
	char path[256];
	sprintf(path, "%s.%d.off", prefix.c_str(), id);
	write_polyhedron(tdbase::make_cube(box), path);
}

void write_box(aab box, const char *path){
	write_polyhedron(tdbase::make_cube(box), path);
}

void write_points(vector<Point> &skeleton, const char *path){
	FILE * fp;
	fp = fopen (path, "w");
	if(fp == NULL){
		log("cannot open %s", path);
		exit(0);
	}

	fprintf(fp, "OFF\n");
	fprintf(fp, "%ld 0 0\n", skeleton.size());
	for(Point &p:skeleton){
		fprintf(fp, "%f %f %f\n",p[0],p[1],p[2]);
	}
	fclose(fp);
}

void write_triangles(vector<Triangle> &triangles, const char *path){
	FILE * fp;
	fp = fopen (path, "w");

	fprintf(fp, "OFF\n");
	fprintf(fp, "%ld %ld 0\n", triangles.size()*3, triangles.size());
	for(Triangle &t:triangles){
		for(int i=0;i<3;i++){
			fprintf(fp, "%f %f %f\n",t[i][0],t[i][1],t[i][2]);
		}
	}
	for(int i=0;i<triangles.size();i++){
		fprintf(fp, "3 %d %d %d\n", i*3, i*3+1,i*3+2);
	}
	fclose(fp);
}

void write_triangles(vector<Triangle *> &triangles, const char *path){
	FILE * fp;
	fp = fopen (path, "w");

	fprintf(fp, "OFF\n");
	fprintf(fp, "%ld %ld 0\n", triangles.size()*3, triangles.size());
	for(Triangle *t:triangles){
		for(int i=0;i<3;i++){
			fprintf(fp, "%f %f %f\n",(*t)[i][0],(*t)[i][1],(*t)[i][2]);
		}
	}
	for(int i=0;i<triangles.size();i++){
		fprintf(fp, "3 %d %d %d\n", i*3, i*3+1,i*3+2);
	}
	fclose(fp);
}

void write_voxels(vector<Voxel *> boxes, const char *path){
	write_polyhedron(make_cubes(boxes), path);
}

void write_polyhedron(Polyhedron *mesh, int id){
	char path[256];
	sprintf(path, "/gisdata/%d.off", id);
	write_polyhedron(mesh, path);
}


/*
 * read polyhedrons or HiMesh
 * */

string read_off_stdin(){
	string input_line;
	getline(std::cin, input_line);
	string whole_mesh = input_line;
	if(input_line.find("|") == std::string::npos){
		whole_mesh += "\n";
		while(getline(std::cin, input_line)){
			whole_mesh += input_line+"\n";
		}
	}else{
		boost::replace_all(whole_mesh, "|", "\n");
	}
	return whole_mesh;
}

HiMesh *read_mesh(bool complete_compression){
	string mesh_str = read_off_stdin();
	HiMesh *mesh = new HiMesh(mesh_str, complete_compression);
	return mesh;
}

HiMesh *read_mesh(const char *path, bool complete_compression){
	string mesh_str = read_file(path);
	HiMesh *mesh = new HiMesh(mesh_str, complete_compression);
	return mesh;
}

vector<HiMesh *> read_meshes(const char *path, size_t maxnum){
	string input_line;
	std::ifstream vfile(path);
	vector<Voxel *> vessel_voxels;
	vector<HiMesh *> ret;
	while(std::getline(vfile, input_line)&&ret.size()<maxnum){
		HiMesh *mesh = new HiMesh(input_line, false);
		ret.push_back(mesh);
	}
	vfile.close();
	return ret;
}

HiMesh *poly_to_mesh(Polyhedron *poly){
	stringstream ss;
	ss<<*poly;
	string str = ss.str();
	HiMesh *mesh = new HiMesh(str, true);
	HiMesh *himesh = new HiMesh(mesh);
	delete mesh;
	return himesh;
}

Polyhedron *parse_polyhedron(string &input){
	stringstream ss;
	ss<<input;
	Polyhedron *poly = new Polyhedron();
	ss>>*poly;
	return poly;
}

Polyhedron *read_polyhedron(const char *path){
	string input = read_file(path);
	return parse_polyhedron(input);
}

Polyhedron *read_polyhedron(){
	string input = read_off_stdin();
	return parse_polyhedron(input);
}

vector<Polyhedron *> read_polyhedrons(const char *path, size_t maxnum){
	string input_line;
	std::ifstream vfile(path);
	vector<Voxel *> vessel_voxels;
	vector<Polyhedron *> ret;
	while(std::getline(vfile, input_line)&&ret.size()<maxnum){
		boost::replace_all(input_line, "|", "\n");
		stringstream ss;
		ss<<input_line;
		Polyhedron *p = new Polyhedron();
		ss >> *p;
		ret.push_back(p);
	}
	vfile.close();
	return ret;
}

}



