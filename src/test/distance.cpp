/*
 * Triangle_Compute.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */


#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <boost/program_options.hpp>

#include "../computational_geometry/TriDist.h"
#include "../PPMC/ppmc.h"
#include "../spatial/spatial.h"

using namespace std;
using namespace hispeed;
namespace po = boost::program_options;

int get_segments(MyMesh *mesh, float *set){

	float *cur_S = set;
	size_t size = 0;
	for(MyMesh::Edge_const_iterator eit = mesh->edges_begin(); eit!=mesh->edges_end(); ++eit){
		Point p1 = eit->vertex()->point();
		Point p2 = eit->opposite()->vertex()->point();
		if(p1==p2){
			continue;
		}
		*cur_S = p1.x();
		cur_S++;
		*cur_S = p1.y();
		cur_S++;
		*cur_S = p1.z();
		cur_S++;
		*cur_S = p2.x();
		cur_S++;
		*cur_S = p2.y();
		cur_S++;
		*cur_S = p2.z();
		cur_S++;
		size++;
	}
	return size;
}

int get_triangles(MyMesh *mesh, float *S){
	float *cur_S = S;
	int size = 0;
	for(MyMesh::Face_iterator fit = mesh->facets_begin(); fit!=mesh->facets_end(); ++fit){
		MyMesh::Halfedge_handle heh = fit->halfedge();
		MyMesh::Halfedge_handle hIt = heh;
		int i = 0;
		do {
			MyMesh::Vertex_handle vh = hIt->vertex();
			Point p = vh->point();
			*cur_S = p.x();
			cur_S++;
			*cur_S = p.y();
			cur_S++;
			*cur_S = p.z();
			cur_S++;
			i++;
			hIt = hIt->next();
		} while (hIt != heh&&i<3);
		size++;

	}
	return size;

}



void TriangleDistance(MyMesh *geom1, MyMesh *geom2, int thread_num){
	long size1 = geom1->size_of_facets(), size2 = geom2->size_of_facets();

	printf("%ld X %ld = %ld\n", size1, size2, size1*size2);
	/*
	 * param S: first triangle set size = s1*3*3*sizeof(float)
	 * param T: second triangle set size = s2*3*3*sizeof(float)
	 * param s1: number of triangles in set S
	 * param s2: number of triangles in set T
	 *
	 * */
	float *S = new float[9*size1];
	float *T = new float[9*size2];


	struct timeval start = get_cur_time();
	size1 = get_triangles(geom1, S);
	size2 = get_triangles(geom2, T);
	printf("loading triangles takes %f milliseconds\n", get_time_elapsed(start));

	start = get_cur_time();

	float min_distance = TriDist_batch(S, T, size1, size2, thread_num);

	printf("get distance take %f miliseconds\n",get_time_elapsed(start));
	printf("min distance is: %f\n", min_distance);

	delete S;
	delete T;
}



void SegmentDistance(MyMesh *geom1, MyMesh *geom2, bool use_gpu, int thread_num){
	long size1 = geom1->size_of_halfedges()/2, size2 = geom2->size_of_halfedges()/2;

	/*
	 * param S: first segment set size = s1*2*3*sizeof(float)
	 * param T: second segment set size = s2*2*3*sizeof(float)
	 * param s1: number of triangles in set S
	 * param s2: number of triangles in set T
	 * */
	float *S = new float[6*size1];
	float *T = new float[6*size2];

	struct timeval start = get_cur_time();
	size1 = get_segments(geom1, S);
	size2 = get_segments(geom2, T);
	printf("%ld X %ld = %ld \n", size1, size2, size1*size2);
	printf("loading triangles takes %f milliseconds\n", get_time_elapsed(start));

	start = get_cur_time();
	float min_distance = 0;
	if(use_gpu){
		min_distance = SegDist_batch_gpu(S, T, size1, size2);
	}else{
		min_distance = SegDist_batch(S, T, size1, size2, thread_num);
	}

	printf("get distance take %f miliseconds\n",get_time_elapsed(start));
	printf("min distance is: %f\n", min_distance);

	delete S;
	delete T;
}

int main(int argc, char **argv){
	int compression_rate = 100;
	int thread_num = std::thread::hardware_concurrency();
	bool use_gpu = false;

	try {
		po::options_description desc("Options");
		desc.add_options()
			("help", "this help message")
			("num_threads,n", po::value<int>(&thread_num),
					"thread number for cpu mode, by default use up all available cpus")
			("gpu,g", "use gpu")
			("compression_rate,r", po::value<int>(&compression_rate),
					"compression rate");
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cerr << desc <<	endl;
			return 0;
		}
		if (vm.count("gpu")) {
			use_gpu = true;
		}
		if(thread_num<0){
			thread_num = 1;
		}
		if(thread_num>std::thread::hardware_concurrency()){
			thread_num = std::thread::hardware_concurrency();
		}
	} catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		exit(0);
	}

	MyMesh *compressed_geom1 = hispeed::read_mesh();
	MyMesh *compressed_geom2 = hispeed::read_mesh();
	struct timeval start = get_cur_time();
	MyMesh *geom1 = hispeed::decompress_mesh(compressed_geom1, compression_rate);
	MyMesh *geom2 = hispeed::decompress_mesh(compressed_geom2, compression_rate);
	SegmentDistance(geom1, geom2, use_gpu, thread_num);

	delete geom1;
	delete geom2;
	delete compressed_geom1;
	delete compressed_geom2;
}
