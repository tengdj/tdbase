/*
 * test.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */

#include <stdio.h>
#include <iostream>
#include <type_traits>
#include "util/util.h"
using namespace std;

int main(int argc, char **argv){


	int count = 0;
	for(int i =0;i<1000000;i++){
		if(hispeed::get_rand_sample(atoi(argv[1]))){
			count++;
		}
	}
	cout<<count<<endl;



}



//        int index = 0;
//        char path[256];
//        for(Skeleton_vertex v : CGAL::make_range(vertices(*skeleton))){
//        	mbb box;
//            for(vertex_descriptor vd : (*skeleton)[v].vertices){
//            	auto p = get(CGAL::vertex_point, tmesh, vd);
//            	box.update(p);
//            	if(index<3){
//            		cout<<(*skeleton)[v].point<<" "<<p<<endl;
//            	}
//            }
//            if(index<3){
//            	cout<<box<<endl;
//            }
//			Polyhedron *pbox = hispeed::make_cube(box);
//			sprintf(path,"offs/%d.off",index++);
//			if(index%5==0){
//				hispeed::print_mesh_file(pbox,path);
//			}
//			delete pbox;
//        }

//void get_skeleton_edges(Skeleton &skeleton){
//	for(Skeleton_edge e : CGAL::make_range(edges(skeleton))){
//		const Point& s = skeleton[source(e, skeleton)].point;
//		const Point& t = skeleton[target(e, skeleton)].point;
//		cout << "2 "<<source(e, skeleton)<<" "<< s << " " <<target(e, skeleton)<<" "<< t << "\n";
//	}
//}

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

