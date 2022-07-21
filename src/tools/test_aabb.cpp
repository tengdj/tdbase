/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include <boost/program_options.hpp>
#include <vector>

#include "../join/SpatialJoin.h"
#include "../spatial/himesh.h"
#include "../index/index.h"

using namespace std;
using namespace hispeed;
namespace po = boost::program_options;

int main(int argc, char **argv){
	struct timeval start = get_cur_time();

	int num = 10000;
	if(argc>3){
		num = atoi(argv[3]);
	}

	Tile *tile1 = new Tile(argv[1], num);
	Tile *tile2 = new Tile(argv[2],1);
	tile1->disable_innerpart();
	tile2->disable_innerpart();
	tile1->retrieve_all();
	tile1->advance_all(100);
	logt("load tiles", start);

	char c;
	std::cin >> c;
	start = get_cur_time();
	for(int i=0;i<tile1->num_objects();i++){
		HiMesh *mesh = tile1->get_mesh(i);
		assert(mesh);
		TriangleTree *tree = mesh->get_aabb_tree_triangle();
		tree->build();
		tree->accelerate_distance_queries();
		//log("indexing %d",i);
	}
	logt("indexed %ld objects",start,tile1->num_objects());
	std::cin >> c;
	start = get_cur_time();

	tile2->retrieve_all();
	tile2->advance_all(100);
	HiMesh *nuc = tile2->get_mesh(0);
	vector<Point> vertices;
	nuc->get_vertices(vertices);
	double mdist = DBL_MAX;
	for(int i=0;i<tile1->num_objects();i++){
		HiMesh *mesh = tile1->get_mesh(i);
		assert(mesh);
		TriangleTree *tree = mesh->get_aabb_tree_triangle();
		for(Point &p:vertices){
			FT sqd = tree->squared_distance(p);
			double dist = (double)CGAL::to_double(sqd);
			mdist = min(mdist, dist);
		}
	}
	logt("querying 1 %f", start, mdist);
	std::cin >> c;
	start = get_cur_time();
	for(int i=0;i<tile1->num_objects();i++){
		HiMesh *mesh = tile1->get_mesh(i);
		assert(mesh);
		TriangleTree *tree = mesh->get_aabb_tree_triangle();
		for(Point &p:vertices){
			FT sqd = tree->squared_distance(p);
			double dist = (double)CGAL::to_double(sqd);
			mdist = min(mdist, dist);
		}
	}
	logt("querying 2 %f", start, mdist);
	std::cin >> c;

}
