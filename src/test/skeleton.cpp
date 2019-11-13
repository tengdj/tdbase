/*
 * skeleton.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: teng
 */


#include "../index/skeleton.h"
#include "../PPMC/ppmc.h"
#include "../util/util.h"

using namespace std;
using namespace hispeed;

int main(int argc, char **argv){

	MyMesh *mesh = read_mesh();
	MyMesh *decompressed = decompress_mesh(mesh, 100);

	std::vector<Point> points;
	struct timeval start = get_cur_time();
	Skeleton *skeleton = extract_skeleton(decompressed);
	cout<<"extracting skeleton takes "<<get_time_elapsed(start)<<" ms"<<endl;
	start  = hispeed::get_cur_time();
	get_skeleton_points(*skeleton, points);
	hispeed::get_skeleton_edges(*skeleton);
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



