/*
 * test.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */

#include <stdio.h>
#include <iostream>
using namespace std;

float TriDist(float S[3][3], float T[3][3]){
	return 0.0;
}

int main(int argc, char **argv){

	int *S = new int[9];
	for(int i=0;i<9;i++){
		S[i]=i;
	}

	for(int i=0;i<3;i++){
		int *curs = S+i*3;
		for(int j=0;j<3;j++){
			cout<<curs[j]<<endl;
		}
	}
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

