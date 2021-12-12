/*
 * decomp.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: teng
 */


#include "../PPMC/ppmc.h"
#include "../util/util.h"
#include "../spatial/spatial.h"
#include "../spatial/himesh.h"
#include <algorithm>
#include "zlib.h"

using namespace hispeed;

int box_num = 600;
int main(int argc, char **argv){

	int start_lod = 0;
	int end_lod = 10;
	if(argc>1){
		start_lod = atoi(argv[1]);
		end_lod = start_lod;
	}
	assert(start_lod>=0&&start_lod<=10);
	// Init the random number generator.
	log("start compressing");
	struct timeval starttime = get_cur_time();

	MyMesh *compressed = read_mesh();
	assert(compressed->size_of_border_edges()&&"must be manifold");
	log("%d vertices %d edges %d faces",compressed->size_of_vertices(), compressed->size_of_halfedges()/2, compressed->true_triangle_size());
	compressed->completeOperation();
	logt("compress", starttime);
	log("start decompressing");

	MyMesh *decompressed = hispeed::decompress_mesh(compressed, 100);
	HiMesh *himesh = new HiMesh(decompressed->p_data,decompressed->dataOffset);

	himesh->advance_to(30);
	float *triangles_low = new float[9*himesh->size_of_facets()];
	size_t low_size = himesh->fill_triangles(triangles_low);

	himesh->advance_to(100);
	float *triangles_high = new float[9*himesh->size_of_facets()];
	size_t high_size = himesh->fill_triangles(triangles_high);
	himesh->computeBoundingBox();
	Point down = himesh->bbMin;
	Point up = himesh->bbMax;

	aab box(down.x(),down.y(),down.z(),up.x(),up.y(),up.z());
	Polyhedron *b = ::make_cube(box);
	write_polyhedron(b, "/gisdata/bbox.OFF");

	vector<float *> tissues;

	float x[8];
	float y[8];
	float z[8];
	for(int i=0;i<box_num;i++){
		double start_x = get_rand_double()*(box.max[0]-box.min[0])+box.min[0];
		double start_y = get_rand_double()*(box.max[1]-box.min[1])+box.min[1];
		double start_z = get_rand_double()*(box.max[2]-box.min[2])+box.min[2];
		double end_x = start_x + (box.max[0]-box.min[0])*(0.01+get_rand_double()/10);
		double end_y = start_y + (box.max[1]-box.min[1])*(0.01+get_rand_double()/10);
		double end_z = start_z + (box.max[2]-box.min[2])*(0.01+get_rand_double()/10);
		x[0] = start_x;
		y[0] = start_y;
		z[0] = start_z;
		x[1] = end_x;
		y[1] = start_y;
		z[1] = start_z;
		x[2] = end_x;
		y[2] = end_y;
		z[2] = start_z;
		x[3] = start_x;
		y[3] = end_y;
		z[3] = start_z;
		x[4] = start_x;
		y[4] = start_y;
		z[4] = end_z;
		x[5] = end_x;
		y[5] = start_y;
		z[5] = end_z;
		x[6] = end_x;
		y[6] = end_y;
		z[6] = end_z;
		x[7] = start_x;
		y[7] = end_y;
		z[7] = end_z;
		float *buf = new float[9*12];
		float *cur = buf;
		*(cur++) = x[0];*(cur++) = y[0];*(cur++) = z[0];
		*(cur++) = x[3]-x[0];*(cur++) = y[3]-y[0];*(cur++) = z[3]-z[0];
		*(cur++) = x[1]-x[0];*(cur++) = y[1]-y[0];*(cur++) = z[1]-z[0];

		*(cur++) = x[1];*(cur++) = y[1];*(cur++) = z[1];
		*(cur++) = x[3]-x[1];*(cur++) = y[3]-y[1];*(cur++) = z[3]-z[1];
		*(cur++) = x[2]-x[1];*(cur++) = y[2]-y[1];*(cur++) = z[2]-z[1];

		*(cur++) = x[0];*(cur++) = y[0];*(cur++) = z[0];
		*(cur++) = x[1]-x[0];*(cur++) = y[1]-y[0];*(cur++) = z[1]-z[0];
		*(cur++) = x[5]-x[0];*(cur++) = y[5]-y[0];*(cur++) = z[5]-z[0];

		*(cur++) = x[0];*(cur++) = y[0];*(cur++) = z[0];
		*(cur++) = x[5]-x[0];*(cur++) = y[5]-y[0];*(cur++) = z[5]-z[0];
		*(cur++) = x[4]-x[0];*(cur++) = y[4]-y[0];*(cur++) = z[4]-z[0];

		*(cur++) = x[1];*(cur++) = y[1];*(cur++) = z[1];
		*(cur++) = x[2]-x[1];*(cur++) = y[2]-y[1];*(cur++) = z[2]-z[1];
		*(cur++) = x[6]-x[1];*(cur++) = y[6]-y[1];*(cur++) = z[6]-z[1];

		*(cur++) = x[1];*(cur++) = y[1];*(cur++) = z[1];
		*(cur++) = x[6]-x[1];*(cur++) = y[6]-y[1];*(cur++) = z[6]-z[1];
		*(cur++) = x[5]-x[1];*(cur++) = y[5]-y[1];*(cur++) = z[5]-z[1];

		*(cur++) = x[2];*(cur++) = y[2];*(cur++) = z[2];
		*(cur++) = x[3]-x[2];*(cur++) = y[3]-y[2];*(cur++) = z[3]-z[2];
		*(cur++) = x[7]-x[2];*(cur++) = y[7]-y[2];*(cur++) = z[7]-z[2];

		*(cur++) = x[2];*(cur++) = y[2];*(cur++) = z[2];
		*(cur++) = x[7]-x[2];*(cur++) = y[7]-y[2];*(cur++) = z[7]-z[2];
		*(cur++) = x[6]-x[2];*(cur++) = y[6]-y[2];*(cur++) = z[6]-z[2];

		*(cur++) = x[5];*(cur++) = y[5];*(cur++) = z[5];
		*(cur++) = x[6]-x[5];*(cur++) = y[6]-y[5];*(cur++) = z[6]-z[5];
		*(cur++) = x[7]-x[5];*(cur++) = y[7]-y[5];*(cur++) = z[7]-z[5];

		*(cur++) = x[4];*(cur++) = y[4];*(cur++) = z[4];
		*(cur++) = x[5]-x[4];*(cur++) = y[5]-y[4];*(cur++) = z[5]-z[4];
		*(cur++) = x[7]-x[4];*(cur++) = y[7]-y[4];*(cur++) = z[7]-z[4];

		*(cur++) = x[0];*(cur++) = y[0];*(cur++) = z[0];
		*(cur++) = x[4]-x[0];*(cur++) = y[4]-y[0];*(cur++) = z[4]-z[0];
		*(cur++) = x[3]-x[0];*(cur++) = y[3]-y[0];*(cur++) = z[3]-z[0];

		*(cur++) = x[3];*(cur++) = y[3];*(cur++) = z[3];
		*(cur++) = x[4]-x[3];*(cur++) = y[4]-y[3];*(cur++) = z[4]-z[3];
		*(cur++) = x[7]-x[3];*(cur++) = y[7]-y[3];*(cur++) = z[7]-z[3];
		tissues.push_back(buf);
		//Polyhedron *bx = ::make_cube(aab(start_x,start_y,start_z,end_x,end_y,end_z));
		//write_polyhedron(bx, i);
	}

	logt("generate random boxes",starttime);

	int low_t = 0;
	for(int i=0;i<box_num;i++){
		low_t += TriInt_single(triangles_low, tissues[i], low_size, 12);
	}
	logt("low resolution %ld %d",starttime,low_size,low_t);

	int high_t = 0;
	for(int i=0;i<box_num;i++){
		high_t += TriInt_single(triangles_high, tissues[i], high_size, 12);
	}
	logt("high resolution %ld %d",starttime, high_size,high_t);


}



