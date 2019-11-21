/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */




#include "../storage/tile.h"
#include "../geometry/TriDist.h"
using namespace std;
using namespace hispeed;

int main(int argc, char **argv){
//	Tile *tile = new Tile("test");
//	SpatialJoin * joiner = new SpatialJoin(tile, tile, JT_nearest);
//	joiner->formalize_computing();

	struct timeval start = get_cur_time();

	MyMesh *geom1 = hispeed::read_mesh();
	MyMesh *geom2 = hispeed::read_mesh();

	printf("parsing geometry takes %f milliseconds\n", get_time_elapsed(start, true));
	size_t size1 = geom1->size_of_halfedges()/2;
	size_t size2 = geom2->size_of_halfedges()/2;
	float *data1 = geom1->get_segments(size1);
	float *data2 = geom2->get_segments(size2);
	printf("%ld X %ld = %ld \n", size1, size2, size1*size2);
	printf("loading segments takes %f milliseconds\n", get_time_elapsed(start, true));

	int *num2 = new int[400];
	int totalnum = 0;
	for(int i=0;i<400;i++){
		num2[i] = hispeed::get_rand_number(4)+1;
		totalnum += num2[i];
	}
	delete num2;
	float *wdata1 = new float[6*400*totalnum];
	float *wdata2 = new float[6*400*totalnum];
	float *cur1 = wdata1;
	float *cur2 = wdata2;
	for(int i=0;i<totalnum;i++){
		memcpy((char *)cur1, (char *)data1, 400*6*sizeof(float));
		memcpy((char *)cur2, (char *)data2, 400*6*sizeof(float));
		cur1 += 400*6;
		cur2 += 400*6;
	}
	printf("copy data takes %f milliseconds\n", get_time_elapsed(start, true));
	float min_distance;
	for(int i=0;i<totalnum;i++){
		min_distance = SegDist_single(wdata1+400*6*i, wdata2+400*6*i, 400, 400);
		//cout<<i<<endl;
	}
	printf("get distance take %f miliseconds\n",get_time_elapsed(start));
	printf("min distance is: %f\n", min_distance);

	delete geom1;
	delete geom2;
	delete data1;
	delete data2;
	delete wdata1;
	delete wdata2;
}
