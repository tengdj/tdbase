/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include <boost/program_options.hpp>
#include <vector>

#include "../join/SpatialJoin.h"
#include "../storage/tile.h"
#include "../spatial/himesh.h"
#include "../index/index.h"

using namespace std;
using namespace hispeed;
namespace po = boost::program_options;

int main(int argc, char **argv){
	struct timeval start = get_cur_time();

	string tile1_path("nuclei_tmp.dt");
	string tile2_path("nuclei_tmp.dt");
	bool use_gpu = false;
	bool intersect = false;
	int num_threads = hispeed::get_num_threads();

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("gpu,g", "compute with GPU")
		("intersect,i", "do intersection instead of join")
		("tile1", po::value<string>(&tile1_path), "path to tile 1")
		("tile2", po::value<string>(&tile2_path), "path to tile 2")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	if(vm.count("gpu")){
		use_gpu = true;
	}
	if(vm.count("intersect")){
		intersect = true;
	}

	Tile *tile1 = new Tile(tile1_path.c_str());
	Tile *tile2 = tile1;
	if(vm.count("tile2")&&tile1_path!=tile2_path){
		tile2 = new Tile(tile2_path.c_str());
	}
	report_time("load tiles", start);

	OctreeNode *tree = tile2->build_octree(2000);
	report_time("build octree", start);



	return 0;


	if(true){
		SpatialJoin *joiner = new SpatialJoin(tile1, tile2);
		if(intersect){
			joiner->intersect(use_gpu, num_threads);
		}else{
			joiner->nearest_neighbor(use_gpu, num_threads);
		}
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

	delete tile1;
	if(tile2!=tile1){
		delete tile2;
	}
}
