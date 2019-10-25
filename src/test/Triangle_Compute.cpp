/*
 * Triangle_Compute.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */


#include <stdio.h>
#include <iostream>
#include <stdlib.h>

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

float get_distance(MyTriangle &t1, MyTriangle &t2, int counter[5]){
	float p[3], q[3], s[3][3], t[3][3];
	for(int i=0;i<3;i++){
		s[i][0] = t1.points[i].x();
		s[i][1] = t1.points[i].y();
		s[i][2] = t1.points[i].z();
		t[i][0] = t2.points[i].x();
		t[i][1] = t2.points[i].y();
		t[i][2] = t2.points[i].z();
	}
	return TriDist(p,q,s,t,counter);
}

int main(int argc, char **argv){

	if(argc>1){
		compression_rate = atoi(argv[1]);
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

	float p[3], q[3], s[3][3], t[3][3];
	long X = 0;
	long Y = 0;
	float min_distance = DBL_MAX-1;

	long size1 = geom1->size_of_facets(), size2 = geom2->size_of_facets();

	printf("%ld X %ld = %ld takes %f milliseconds\n", size1, size2, size1*size2, get_time_elapsed(start));
	start = get_cur_time();
	MyTriangle *tset1 = new MyTriangle[size1];
	MyTriangle *tset2 = new MyTriangle[size2];

	int index = 0;
	for(MyMesh::Face_iterator fit = geom1->facets_begin(); fit!=geom1->facets_end(); ++fit){

		MyMesh::Halfedge_handle heh = fit->halfedge();
		MyMesh::Halfedge_handle hIt = heh;
		int i = 0;
		do {
			MyMesh::Vertex_handle vh = hIt->vertex();
			tset1[index].points[i] = vh->point();
			i++;
			hIt = hIt->next();
		} while (hIt != heh&&i<3);
		index++;
	}

	index = 0;
	for(MyMesh::Face_iterator fit = geom2->facets_begin(); fit!=geom2->facets_end(); ++fit){
		MyMesh::Halfedge_handle heh = fit->halfedge();
		MyMesh::Halfedge_handle hIt = heh;
		int i = 0;
		do {
			MyMesh::Vertex_handle vh = hIt->vertex();
			tset2[index].points[i] = vh->point();
			i++;
			hIt = hIt->next();
		} while (hIt != heh&&i<3);
		index++;
	}
	printf("loading triangles takes %f milliseconds\n", get_time_elapsed(start));

	int counter[5] = {0,0,0,0,0};
	start = get_cur_time();
	for(int i=0;i<size1;i++){
		for(int j=0;j<size2;j++){
			float distance = get_distance(tset1[i],tset2[j],counter);
			if(distance<min_distance){
				min_distance = distance;
			}
			if(distance == 0){
				cout<<i<<" "<<j<<endl;
				for(int k = 0;k<3;k++){
					cout<<tset1[i].points[k]<<"   "<<tset2[j].points[k]<<endl;
				}
				cout<<endl;
			}
		}
	}



	printf("get distance take %f miliseconds\n",get_time_elapsed(start));
	printf("min distance is: %f\n", min_distance);
	printf("stats: %d %d %d %d %d\n",counter[0],counter[1],counter[2],counter[3],counter[4]);


	delete geom1;
	delete geom2;
	delete compressed_geom1;
	delete compressed_geom2;
	delete tset1;
	delete tset2;

}
