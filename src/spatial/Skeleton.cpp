/*
 * Skeleton.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: teng
 */

#include "spatial.h"


namespace hispeed{

void get_skeleton_points(Skeleton &skeleton, std::vector<Point> &P){
    BOOST_FOREACH(Skeleton_vertex v, boost::vertices(skeleton)){
		auto p = skeleton[v].point;
		for(vertex_descriptor vd : skeleton[v].vertices){
			uint idx = vd;
			cout <<v<< " " << idx << "\n";
		}
		P.push_back(Point(p.x(),p.y(),p.z()));
    }
}

void get_skeleton_edges(Skeleton &skeleton){
	for(Skeleton_edge e : CGAL::make_range(edges(skeleton))){
		const Point& s = skeleton[source(e, skeleton)].point;
		const Point& t = skeleton[target(e, skeleton)].point;
		cout << "2 "<<source(e, skeleton)<<" "<< s << " " <<target(e, skeleton)<<" "<< t << "\n";
	}
}

Skeleton * extract_skeleton(MyMesh *currentMesh){

#ifdef DEBUG
        cerr << "extracting the Skeleton (advanced)!" << endl;
#endif
        std::stringstream os;
        os << *currentMesh;
        Triangle_mesh tmesh;
        os >> tmesh;
        Skeleton *skeleton = new Skeleton();
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
#ifdef DEBUG
        std::cerr<<"extracted"<<endl;
#endif

        int index = 0;
        char path[256];
        for(Skeleton_vertex v : CGAL::make_range(vertices(*skeleton))){
        	mbb box;
            for(vertex_descriptor vd : (*skeleton)[v].vertices){
            	auto p = get(CGAL::vertex_point, tmesh, vd);
            	box.update(p);
            }
			Polyhedron *pbox = hispeed::make_cube(box);
			sprintf(path,"offs/%d.off",index++);
			if(index%5==0){
				hispeed::print_mesh_file(pbox,path);
			}
			delete pbox;
        }
//    	BOOST_FOREACH(Skeleton_vertex v, boost::vertices(*skeleton)){
//
//    		for(vertex_descriptor vd : (*skeleton)[v].vertices){
//    			box.update(points[vd]);
//    		}

//    	}
        return skeleton;
}


}
