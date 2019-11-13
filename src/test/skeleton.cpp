/*
 * skeleton.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: teng
 */


#include "../spatial/spatial.h"
#include "../PPMC/ppmc.h"
#include "../util/util.h"

using namespace std;
using namespace hispeed;
using namespace CGAL;

int main(int argc, char **argv){

	MyMesh *mesh = read_mesh();
	MyMesh *decompressed = decompress_mesh(mesh, 100);

	std::vector<Point> points;
	struct timeval start = get_cur_time();
	Skeleton *skeleton = extract_skeleton(decompressed);
	cout<<"extracting skeleton takes "<<get_time_elapsed(start)<<" ms"<<endl;
	start  = hispeed::get_cur_time();
//	get_skeleton_points(*skeleton, points);
//	auto p = CGAL::get(CGAL::vertex_point, decompressed, 1);
	for(Vertex_iterator v = decompressed->vertices_begin(); v != decompressed->vertices_end(); ++v){
		points.push_back(v->point());
	}

//	int index = 0;
//	char path[256];
//	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
//
//		mbb box;
//		for(vertex_descriptor vd : (*skeleton)[v].vertices){
//			box.update(points[vd]);
//		}
//		Polyhedron *pbox = hispeed::make_cube(box);
//		std::stringstream os;
//		sprintf(path,"offs/%d.off",index++);
//		if(index%3==0){
//			hispeed::print_mesh_file(pbox,path);
//		}
//		delete pbox;
//	}


	//hispeed::get_skeleton_edges(*skeleton);
	cout<<"get skeleton points takes "<<get_time_elapsed(start)<<" ms"<<endl;

//	int i=0;
//	for(Point p: points){
//		cout<<p.x()<<endl;
//	}
//	cout<<endl;
//	for(Point p: points){
//		cout<<p.y()<<endl;
//	}
//	cout<<endl;
//	for(Point p: points){
//		cout<<p.z()<<endl;
//	}
//	cout<<endl;
	points.clear();

	delete decompressed;
	delete mesh;
	delete skeleton;

}



