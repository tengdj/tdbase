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

#include "../triangle/TriDist.h"
#include "spatial.h"
#include "../PPMC/mymesh.h"
#include "../PPMC/configuration.h"
#include "PPMC/ppmc.h"

using namespace std;

int compression_rate = 100;
MyMesh *get_mesh(string input_line){
	int i_mode = COMPRESSION_MODE_ID; // compression mode
	unsigned i_quantBit = 12;
	unsigned i_decompPercentage = 100;
	bool b_allowConcaveFaces = true;
	// Init the random number generator.
	srand(PPMC_RANDOM_CONSTANT);
	MyMesh *compressed = new MyMesh(i_decompPercentage,
				 i_mode, i_quantBit,
				 b_allowConcaveFaces,
				 input_line.c_str(), input_line.size());
	compressed->completeOperation();
	return compressed;
}

MyMesh *decompress(MyMesh *compressed){
	int i_mode = DECOMPRESSION_MODE_ID; // compression mode
	unsigned i_quantBit = 12;
	unsigned i_decompPercentage = compression_rate;
	bool b_allowConcaveFaces = true;

	// Init the random number generator.
	srand(PPMC_RANDOM_CONSTANT);
	MyMesh *decompressed = new MyMesh(i_decompPercentage,
				 i_mode, i_quantBit,
				 b_allowConcaveFaces,
				 compressed->p_data, compressed->dataOffset);
	decompressed->completeOperation();
	return decompressed;
}

typedef struct MyTriangle {
	Point points[3];
} MyTriangle;

//float get_distance(MyTriangle &t1, MyTriangle &t2){
//	float p[3], q[3], s[3][3], t[3][3];
//	for(int i=0;i<3;i++){
//		s[i][0] = t1.points[i].x();
//		s[i][1] = t1.points[i].y();
//		s[i][2] = t1.points[i].z();
//		t[i][0] = t2.points[i].x();
//		t[i][1] = t2.points[i].y();
//		t[i][2] = t2.points[i].z();
//	}
//	return TriDist(p,q,s,t);
//}

int main(int argc, char **argv){

	int thread_num = std::thread::hardware_concurrency();
	if(argc>1){
		thread_num = atoi(argv[1]);
	}
	string input_line1, input_line2;
	getline(std::cin, input_line1);
	boost::replace_all(input_line1, "|", "\n");
	getline(std::cin, input_line2);
	boost::replace_all(input_line2, "|", "\n");

	MyMesh *compressed_geom1 = get_mesh(input_line1);
	MyMesh *compressed_geom2 = get_mesh(input_line2);
	struct timeval start = get_cur_time();
	MyMesh *geom1 = decompress(compressed_geom1);
	MyMesh *geom2 = decompress(compressed_geom2);


	long size1 = geom1->size_of_facets(), size2 = geom2->size_of_facets();

	printf("%ld X %ld = %ld takes %f milliseconds\n", size1, size2, size1*size2, get_time_elapsed(start));
	/*
	 * param dist: head address of space for distances, size = s1*s2*sizeof(float)
	 * param P: closest points in S size = 3*s1*s2
	 * param Q: closest points in T size = 3*s1*s2
	 * param S: first triangle set size = s1*3*3*sizeof(float)
	 * param T: second triangle set size = s2*3*3*sizeof(float)
	 * param s1: number of triangles in set S
	 * param s2: number of triangles in set T
	 *
	 * */
	float *S = new float[9*size1];
	float *T = new float[9*size2];


	start = get_cur_time();

	float *cur_S = S;
	for(MyMesh::Face_iterator fit = geom1->facets_begin(); fit!=geom1->facets_end(); ++fit){

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
	}

	float *cur_T = T;
	for(MyMesh::Face_iterator fit = geom2->facets_begin(); fit!=geom2->facets_end(); ++fit){
		MyMesh::Halfedge_handle heh = fit->halfedge();
		MyMesh::Halfedge_handle hIt = heh;
		int i = 0;
		do {
			MyMesh::Vertex_handle vh = hIt->vertex();
			Point p = vh->point();
			*cur_T = p.x();
			cur_T++;
			*cur_T = p.y();
			cur_T++;
			*cur_T = p.z();
			cur_T++;
			i++;
			hIt = hIt->next();
		} while (hIt != heh&&i<3);
	}
	printf("loading triangles takes %f milliseconds\n", get_time_elapsed(start));

	start = get_cur_time();

	float min_distance = TriDist_batch(S, T, size1, size2, thread_num);


	printf("get distance take %f miliseconds\n",get_time_elapsed(start));
	printf("min distance is: %f\n", min_distance);


	delete geom1;
	delete geom2;
	delete compressed_geom1;
	delete compressed_geom2;

	delete S;
	delete T;

}
