/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */



#include "../join/SpatialJoin.h"
#include "../storage/tile.h"
#include "../spatial/himesh.h"
#include <vector>
using namespace std;
using namespace hispeed;

int main(int argc, char **argv){
	struct timeval start = get_cur_time();

	Tile *tile1 = new Tile("nuclei.dt");
	Tile *tile2 = new Tile("vessel.dt");
	report_time("load tiles", start);
	if(true){
		SpatialJoin *joiner = new SpatialJoin(tile1, tile1);
		//joiner->nearest_neighbor(true, 1);
		joiner->intersect(false);
		report_time("total join", start);
		delete joiner;
	}

	if(false){
		HiMesh *mesh = tile2->get_mesh(0,100);
		HiMesh *mesh2 = tile2->get_mesh(1,100);
		report_time("get meshs", start);
		SegTree *tree = mesh->get_aabb_tree();
		report_time("get aabb tree", start);
		std::vector<Point> points;
		mesh2->get_vertices(points);
		float mindist = DBL_MAX;
		for(Point p:points){
			float dist = (float)CGAL::to_double(tree->squared_distance(p));
			if(dist<mindist){
				mindist = dist;
			}
		}
		cout<<sqrt(mindist)<<endl;
		report_time("getting min", start);
	}

	if(false){
		TriangleTree *tree1 = get_aabb_tree(hispeed::read_polyhedron());
		TriangleTree *tree2 = get_aabb_tree(hispeed::read_polyhedron());
		report_time("get aabb tree", start);
		for(int i=0;i<tile1->num_objects();i++){
			float mindist1 = DBL_MAX;
			float mindist2 = DBL_MAX;
			HiMesh *mesh = tile1->get_mesh(i,100);
			std::vector<Point> points;
			mesh->get_vertices(points);
			for(Point p:points){
				float dist = (float)CGAL::to_double(tree1->squared_distance(p));
				if(dist<mindist1){
					mindist1 = dist;
				}
				dist = (float)CGAL::to_double(tree2->squared_distance(p));
				if(dist<mindist2){
					mindist2 = dist;
				}
			}
			points.clear();
		}
		report_time("getting min", start);
	}


//	report_time("building aabb tree", start);
//	cout<<(float)CGAL::to_double(tree->squared_distance(Point(0, 0, 0)))<<endl;

	delete tile1;
	delete tile2;
}
