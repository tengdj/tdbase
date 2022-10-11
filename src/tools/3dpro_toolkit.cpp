/*
 * 3dpro_tool.cpp
 *
 *  Created on: Mar 29, 2022
 *      Author: teng
 */

#include <algorithm>
#include <zlib.h>

#include "himesh.h"
#include "tile.h"
#include "../PPMC/ppmc.h"
#include "util.h"

using namespace hispeed;

/*
 * print himesh to wkt
 *
 * */
static void himesh_to_wkt(int argc, char **argv){

	Tile *tile = new Tile(argv[1]);
	tile->retrieve_all();
	tile->advance_all(100);
	for(int i=0;i<tile->num_objects();i++){
		cout<<i<<"|"<<tile->get_mesh(i)->to_wkt()<<endl;
	}
	delete tile;
}


/*
 * profiling the distribution of protruding vertices
 *
 * */
static void profile_protruding(int argc, char **argv){

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
	log("%d vertices %d edges %d faces",compressed->size_of_vertices(), compressed->size_of_halfedges()/2, compressed->size_of_facets());
	compressed->profileProtruding();
}

/*
 * get the voxel boxes
 * */
static void get_voxel_boxes(int argc, char **argv){
	struct timeval start = get_cur_time();
	MyMesh *mesh = hispeed::read_mesh();
	mesh->completeOperation();
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);
	int voxel_num = 100;
	if(argc>=2){
		voxel_num = atoi(argv[1]);
	}
	vector<Voxel *> voxels = himesh->generate_voxels_skeleton(voxel_num);

	float vol = 0.0;
	for(Voxel *v:voxels){
		vol += v->volume();
	}

	log("%ld voxels are generated %f volumn", voxels.size(), vol);
	write_voxels(voxels, "/gisdata/skeleton_voxels.off");
	hispeed::write_box(himesh->get_box(), "/gisdata/aab.off");
	delete himesh;
}


/*
 * get the skeleton points
 * */
static void get_skeleton(int argc, char **argv){
	assert(argc>2);
	struct timeval start = get_cur_time();
	MyMesh *mesh = read_mesh(argv[1]);
	mesh->completeOperation();
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);

	int voxel_num = 100;
	if(argc>3){
		voxel_num = atoi(argv[3]);
	}

	vector<Point> skeleton = himesh->get_skeleton_points(voxel_num);
	hispeed::write_points(skeleton, argv[2]);

	log("%ld points in the skeleton", skeleton.size());
	skeleton.clear();
	delete mesh;
	delete himesh;
}

/*
 * get the voxel boxes
 * */
static void voxelize(int argc, char **argv){
	assert(argv>2);
	struct timeval start = get_cur_time();
	MyMesh *mesh = read_mesh(argv[1]);

	mesh->completeOperation();
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);
	int voxel_num = 100;
	if(argc>3){
		voxel_num = atoi(argv[3]);
	}
	vector<Voxel *> voxels = himesh->voxelization(voxel_num);
	float vol = 0.0;
	for(Voxel *v:voxels){
		vol += v->volume();
	}

	log("%ld voxels are generated %f volumn", voxels.size(), vol);

	write_voxels(voxels, argv[2]);
	delete himesh;
}

/*
 * get the distances between the objects in two tiles with different LODs
 *
 * */
static void profile_distance(int argc, char **argv){

	Tile *tile1 = new Tile("/home/teng/git/3DPro/src/teng_v_nv1_nu200_s10_vs100_r30.dt");
	Tile *tile2 = new Tile("/home/teng/git/3DPro/src/teng_n_nv1_nu200_s10_vs100_r30.dt");

	for(int i=10;i<=100;i+=10){
		tile1->advance_all(i);
		tile2->advance_all(i);
		float *triangles_v = new float[9*tile1->get_mesh(0)->size_of_facets()];
		tile1->get_mesh(0)->fill_triangles(triangles_v);
		double dist = 0;
		for(int j=0;j<tile2->num_objects();j++) {
			float *triangles_n = new float[9*tile2->get_mesh(j)->size_of_facets()];
			tile2->get_mesh(j)->fill_triangles(triangles_n);
			size_t s1 = tile1->get_mesh(0)->size_of_facets();
			size_t s2 = tile2->get_mesh(j)->size_of_facets();
			float d = TriDist_single(triangles_v,triangles_n,s1,s2);
			dist += sqrt(d);
			delete []triangles_n;
		}
		delete []triangles_v;
		printf("%d %f\n",i,dist/tile2->num_objects());
	}

	delete tile1;
	delete tile2;
}

/*
 * profiling the performance of decoding
 * */
void profile_decoding(int argc, char **argv){

	int start_lod = 0;
	int end_lod = 10;
	if(argc>2){
		start_lod = atoi(argv[2]);
		end_lod = start_lod;
	}
	assert(start_lod>=0&&start_lod<=10);
	// Init the random number generator.
	log("start compressing");
	MyMesh *compressed = read_mesh(argv[1]);

	struct timeval starttime = get_cur_time();
	//assert(compressed->size_of_border_edges()&&"must be manifold");
	compressed->completeOperation();
	logt("compress", starttime);

	const int itertime = 10;

	MyMesh *testc[itertime];
	for(int i=0;i<itertime;i++){
		testc[i] = read_mesh(argv[1]);
	}
	struct timeval sst = get_cur_time();
	for(int i=0;i<itertime;i++){
		testc[i]->completeOperation();
	}
	log("compress %.4f",get_time_elapsed(sst)/itertime);

	log("%d vertices %d edges %d faces",compressed->size_of_vertices(), compressed->size_of_halfedges()/2, compressed->true_triangle_size());

	log("start decompressing");

	HiMesh *himesh;
	char path[256];
	for(int i=start_lod;i<=end_lod;i++){
		int lod = 10*i;
		MyMesh *tested[itertime];
		for(int t=0;t<itertime;t++){
			tested[t] = hispeed::decompress_mesh(compressed, lod, false);
		}
		starttime = get_cur_time();
		for(int t=0;t<itertime;t++){
			tested[t]->completeOperation();
		}
		double testedtime = get_time_elapsed(starttime,true);

		for(int t=0;t<itertime;t++){
			delete tested[t];
		}
		MyMesh *decompressed = hispeed::decompress_mesh(compressed, lod);
		decompressed->completeOperation();
		logt("decompress %3d lod %5d vertices %5d edges %5d faces avg(%.4f)", starttime, lod,
				decompressed->size_of_vertices(), decompressed->size_of_halfedges()/2, decompressed->true_triangle_size(), testedtime/itertime);
		sprintf(path,"/gisdata/lod.%d.off", lod);
		decompressed->writeMeshOff(path);
		if(lod==100){
			himesh = new HiMesh(decompressed->p_data,decompressed->dataOffset);
		}
		delete decompressed;
	}
	delete compressed;
	starttime = get_cur_time();
	float *vertices = NULL;
	himesh->advance_to(100);
	logt("decompress", starttime);
	size_t size = himesh->fill_vertices(vertices);
	logt("fill vertices %d with %d bytes (%ld bytes)", starttime, size, size*3*sizeof(float),himesh->dataOffset);
	char *zcomp = new char[size*3*sizeof(float)];
	unsigned long compressedsize;
	for(int i=1;i<10;i++){
		int nResult = compress2((unsigned char *)zcomp, &compressedsize, (unsigned char *)vertices, size*3*sizeof(float),i);
		logt("compress %d level %ld bytes",starttime,i,compressedsize);
	}

	aab box = himesh->get_box();
	hispeed::write_box(box, "/gisdata/box.off");

	delete himesh;
	delete []vertices;

}

/*
 * adjust the size and the position of a polyhedron
 * */
static void adjust_polyhedron(int argc, char **argv){
	if(argc<=3){
		cout<<"usage: adjust shift shrink output_path"<<endl;
		return;
	}
	int sft = atoi(argv[1]);
	int shift[3] = {sft,sft,sft};
	float shrink = atoi(argv[2]);;
	Polyhedron *poly = hispeed::read_polyhedron();
	Polyhedron apoly = adjust_polyhedron(shift,shrink,poly);
	hispeed::write_polyhedron(&apoly, argv[3]);
	delete poly;
}

static void triangulate(int argc, char **argv){
	MyMesh *mesh = read_mesh(argv[1]);
	mesh->completeOperation();
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);

	Polyhedron *poly = himesh->to_triangulated_polyhedron();
	hispeed::write_polyhedron(poly, argv[2]);

	delete poly;
	delete mesh;
	delete himesh;
}

static void compress(int argc, char **argv){
	MyMesh *mesh = read_mesh(argv[1]);
	mesh->writeMeshOff("/gisdata/origin.off");
	assert(mesh);
	mesh->completeOperation();

	mesh->writeMeshOff("/gisdata/compressed.off");

	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	int lod = 100;
	if(argc>2){
		lod = atoi(argv[2]);
	}
	for(int i=0;i<=lod;i+=10){
		himesh->advance_to(lod);
		himesh->writeCurrentOperationMesh("/gisdata/compressed", i);
	}
	delete mesh;
	delete himesh;
}

static void test(int argc, char **argv){

	float point[3] = {0,0,0};
	float triangle[9] = {1,0,0,0,1,0,2,0,0};
	log("%f", PointTriangleDist(point, triangle));

}

int main(int argc, char **argv){
	global_ctx = parse_args(argc, argv);
	if(argc==1){
		cout<<"usage: 3dpro function [args]"<<endl;
		exit(0);
	}
	if(strcmp(argv[1],"himesh_to_wkt") == 0){
		himesh_to_wkt(argc-1,argv+1);
	}else if(strcmp(argv[1],"profile_protruding") == 0){
		profile_protruding(argc-1,argv+1);
	}else if(strcmp(argv[1],"get_voxel_boxes") == 0){
		get_voxel_boxes(argc-1,argv+1);
	}else if(strcmp(argv[1],"profile_distance") == 0){
		profile_distance(argc-1,argv+1);
	}else if(strcmp(argv[1],"profile_decoding") == 0){
		profile_decoding(argc-1,argv+1);
	}else if(strcmp(argv[1],"adjust_polyhedron") == 0){
		adjust_polyhedron(argc-1,argv+1);
	}else if(strcmp(argv[1],"skeleton") == 0){
		get_skeleton(argc-1,argv+1);
	}else if(strcmp(argv[1],"voxelize") == 0){
		voxelize(argc-1,argv+1);
	}else if(strcmp(argv[1],"test") == 0){
		test(argc-1,argv+1);
	}else if(strcmp(argv[1],"compress") == 0){
		compress(argc-1,argv+1);
	}else if(strcmp(argv[1],"triangulate") == 0){
		triangulate(argc-1,argv+1);
	}else{
		cout<<"usage: 3dpro himesh_to_wkt|profile_protruding|get_voxel_boxes|profile_distance|profile_decoding|adjust_polyhedron|skeleton|voxelize [args]"<<endl;
		exit(0);
	}
	return 0;
}

