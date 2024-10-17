/*
 * 3dpro_tool.cpp
 *
 *  Created on: Mar 29, 2022
 *      Author: teng
 */

#include <algorithm>
#include <zlib.h>
#include <thread>

#include "SpatialJoin.h"
#include "himesh.h"
#include "tile.h"
#include "util.h"

using namespace std;
using namespace tdbase;

namespace tdbase{
/*
 * print himesh to wkt
 *
 * */
static void to_wkt(int argc, char **argv){

	Tile *tile = new Tile(argv[1]);
	const char *prefix = argv[2];
	tile->decode_all(100);
	char path[256];
	for(int i=0;i<tile->num_objects();i++){
		sprintf(path, "%s_%d.OFF", prefix, i);
		tile->get_mesh(i)->write_to_off(path);
	}
	delete tile;
}

/*
 * print himesh to sql
 *
 * */
static void to_sql(int argc, char **argv){

	Tile *tile = new Tile(argv[1]);
	char table[100];
	for(int i=20;i<=100;i+=20){
		tile->decode_all(i);
		char path[256];
		sprintf(table, "%s_%d",argv[2],i);
		sprintf(path, "%s.sql",table);
		tile->dump_sql(path, table);
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
	HiMesh *compressed = read_mesh();
	assert(compressed->size_of_border_edges()&&"must be manifold");
	log("%d vertices %d edges %d faces",compressed->size_of_vertices(), compressed->size_of_halfedges()/2, compressed);
	compressed->profileProtruding();
}

/*
 * get the voxel boxes
 * */
static void get_voxel_boxes(int argc, char **argv){
	struct timeval start = get_cur_time();
	HiMesh *mesh = read_mesh(argv[1], false);
	int voxel_size = 100;
	if(argc>=3){
		voxel_size = atoi(argv[2]);
	}
	vector<Voxel *> voxels = mesh->generate_voxels_skeleton(voxel_size);

	float vol = 0.0;
	for(Voxel *v:voxels){
		vol += v->volume();
	}

	logt("%ld voxels are generated %f volumn", start, voxels.size(), vol);

	write_voxels(voxels, "/share/skeleton_voxels.off");
	tdbase::write_box(mesh->get_mbb(), "/share/aab.off");
	delete mesh;
}


/*
 * get the skeleton points
 * */
static void get_skeleton(int argc, char **argv){
	assert(argc>2);
	struct timeval start = get_cur_time();
	HiMesh *mesh = read_mesh(argv[1], false);
	int voxel_num = 100;
	if(argc>3){
		voxel_num = atoi(argv[3]);
	}
	logt("load mesh", start);
	vector<Point> skeleton = mesh->get_skeleton_points(voxel_num);
	logt("get skeleton", start);

	tdbase::write_points(skeleton, argv[2]);

	log("%ld points in the skeleton", skeleton.size());
	skeleton.clear();
	delete mesh;
}

/*
 * get the voxel boxes
 * */
static void voxelize(int argc, char **argv){
	assert(argc>2);
	struct timeval start = get_cur_time();
	HiMesh *mesh = read_mesh(argv[1], false);
	int voxel_num = 100;
	if(argc>3){
		voxel_num = atoi(argv[3]);
	}
	vector<Voxel *> voxels = mesh->voxelization(voxel_num);
	float vol = 0.0;
	for(Voxel *v:voxels){
		vol += v->volume();
	}

	log("%ld voxels are generated %f volumn", voxels.size(), vol);

	write_voxels(voxels, argv[2]);
	delete mesh;
}

/*
 * get the distances between the objects in two tiles with different LODs
 *
 * */
static void profile_distance(int argc, char **argv){

	Tile *tile1 = new Tile(argv[1]);
	Tile *tile2 = new Tile(argv[2]);

	for(int i=10;i<=100;i+=10){
		tile1->decode_all(i);
		tile2->decode_all(i);
		double dist = 0;
		for(int j=0;j<tile2->num_objects();j++) {
			float d = tile1->get_mesh(0)->distance(tile2->get_mesh(j));
			dist += sqrt(d);
		}
		printf("%d %f\n",i,dist/tile2->num_objects());
	}

	delete tile1;
	delete tile2;
}

/*
 * profiling the performance of decoding
 * */
static void profile_decoding(int argc, char **argv){

	int start_lod = 0;
	int end_lod = 10;
	if(argc>2){
		start_lod = atoi(argv[2]);
		end_lod = start_lod;
	}
	assert(start_lod>=0&&start_lod<=10);
	// Init the random number generator.
	log("start compressing");
	struct timeval starttime = get_cur_time();
	HiMesh *compressed = read_mesh(argv[1], true);
	//assert(compressed->size_of_border_edges()&&"must be manifold");
	logt("compress", starttime);

	const int itertime = 5;

	HiMesh *testc[itertime];
	struct timeval sst = get_cur_time();
	for(int i=0;i<itertime;i++){
		testc[i] = read_mesh(argv[1], true);
	}
	logt("compress to %d vertices %d edges %d faces",sst,compressed->size_of_vertices(), compressed->size_of_halfedges()/2, compressed->size_of_triangles());

	log("start decompressing");

	char path[256];
	for(int i=start_lod;i<=end_lod;i++){
		int lod = 10*i;
		HiMesh *tested[itertime];
		for(int t=0;t<itertime;t++){
			tested[t] = new HiMesh(compressed);
		}
		starttime = get_cur_time();
		for(int t=0;t<itertime;t++){
			tested[t]->decode(lod);
		}
		double testedtime = get_time_elapsed(starttime,true);
		logt("decompress %3d lod %5d vertices %5d edges %5d faces avg(%.4f)", starttime, lod,
				tested[0]->size_of_vertices(), tested[0]->size_of_halfedges()/2, tested[0]->size_of_triangles(), testedtime/itertime);
		sprintf(path,"/gisdata/lod.%d.off", lod);
		tested[0]->write_to_off(path);
		for(int t=0;t<itertime;t++){
			delete tested[t];
		}
	}
	starttime = get_cur_time();
	HiMesh *himesh = new HiMesh(compressed);
	himesh->decode(100);
	logt("decompress", starttime);
	float *vertices = NULL;
	size_t size = himesh->fill_vertices(vertices);
	logt("fill vertices %d with %d bytes (%ld bytes)", starttime, size, size*3*sizeof(float),himesh->get_data_size());
	char *zcomp = new char[size*3*sizeof(float)];
	unsigned long compressedsize;
	for(int i=1;i<10;i++){
		int nResult = compress2((unsigned char *)zcomp, &compressedsize, (unsigned char *)vertices, size*3*sizeof(float),i);
		logt("compress %d level %ld bytes",starttime,i,compressedsize);
	}

	aab box = himesh->get_mbb();
	tdbase::write_box(box, "/gisdata/box.off");

	delete himesh;
	delete compressed;
	delete []vertices;
}

static void aabb(int argc, char **argv){
	struct timeval start = get_cur_time();

	int num = 10000;
	if(argc>3){
		num = atoi(argv[3]);
	}

	Tile *tile1 = new Tile(argv[1], num);
	Tile *tile2 = new Tile(argv[2],1);
	tile1->decode_all(100);
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

	tile2->decode_all(100);
	HiMesh *nuc = tile2->get_mesh(0);
	list<Point> vertices = nuc->get_vertices();
	double mdist = DBL_MAX;
	for(int i=0;i<tile1->num_objects();i++){
		HiMesh *mesh = tile1->get_mesh(i);
		assert(mesh);
		TriangleTree *tree = mesh->get_aabb_tree_triangle();
		for(Point &p:vertices){
			FT sqd = tree->squared_distance(p);
			double dist = (double)CGAL::to_double(sqd);
			mdist = std::min(mdist, dist);
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
			mdist = std::min(mdist, dist);
		}
	}
	logt("querying 2 %f", start, mdist);
	std::cin >> c;

}


/*
 * adjust the size and the position of a polyhedron
 * */
static void adjust_polyhedron(int argc, char **argv){
	if(argc<=3){
		cout<<"usage: adjust shift shrink output_path"<<endl;
		return;
	}
	float sft = atof(argv[1]);
	float shrink = atof(argv[2]);
	HiMesh *mesh = tdbase::read_mesh();
	char path[256];
	sprintf(path, "%s/original.off", argv[3]);
	mesh->write_to_off(path);
	mesh->shift(sft, sft, sft);
	mesh->shrink(shrink);
	sprintf(path, "%s/adjusted_%.2f_%.2f.off",argv[3],sft,shrink);
	mesh->write_to_off(path);
	delete mesh;
}

static void triangulate(int argc, char **argv){
	HiMesh *mesh = read_mesh(argv[1], false);
	Polyhedron *poly = mesh->to_triangulated_polyhedron();
	tdbase::write_polyhedron(poly, argv[2]);

	delete poly;
	delete mesh;
}

static void sample(int argc, char **argv){

	HiMesh::sampling_rate = atoi(argv[3]);
	log("sampling rate %d", HiMesh::sampling_rate);

	struct timeval start = get_cur_time();
	HiMesh *mesh = read_mesh(argv[1]);

	unordered_set<Point> point_set;
	mesh->sample_points(point_set);
	vector<Point> points;
	points.assign(point_set.begin(), point_set.end());
	tdbase::write_points(points, argv[2]);
	delete mesh;
}

static void simplify(int argc, char **argv){


}

static void compress(int argc, char **argv){

	std::string str(argv[1]);

	global_ctx.verbose = 2;
	if(argc>2){
		HiMesh::sampling_rate = atoi(argv[2]);
		log("%d",HiMesh::sampling_rate);
	}

	if(argc>3){
		int cm = atoi(argv[3]);
		HiMesh::calculate_method = (tdbase::Hausdorff_Computing_Type)cm;
	}
	struct timeval start = get_cur_time();
	HiMesh *mesh = read_mesh(argv[1], true);
	logt("compress", start);
	HiMesh *hm = new HiMesh(mesh);
	int lod = 100;

	char path[256];
	std::string base_filename = str.substr(str.find_last_of("/\\") + 1, str.find_last_of('.') -str.find_last_of("/\\") - 1);

	for(uint32_t i=20;i<=lod;i+=20){
		hm->decode(i);
		logt("decode to %d with %d vertices", start, i, hm->size_of_vertices());
		//log("%d %f", i, himesh->getHausdorfDistance());
	    sprintf(path, "%s_%d.off", base_filename.c_str(), i);
	    hm->write_to_off(path);
	}

	delete mesh;
	delete hm;
}

static void distance(int argc, char **argv){
	assert(argc>4);
	int n1 = atoi(argv[3]);
	int n2 = atoi(argv[4]);

	Tile *tile1 = new Tile(argv[1],n1+1);
	Tile *tile2 = new Tile(argv[2],n2+1);

	int lod = 100;
	if(argc>5){
		lod = atoi(argv[5]);
	}

	tile1->decode_all(lod);
	tile2->decode_all(lod);

	HiMesh *mesh1 = tile1->get_mesh(n1);
	HiMesh *mesh2 = tile2->get_mesh(n2);
	log("%f %f", mesh1->distance(mesh2), mesh1->distance_tree(mesh2));

	delete tile1;
	delete tile2;
}

static void intersect(int argc, char **argv){
	assert(argc>4);
	int n1 = atoi(argv[3]);
	int n2 = atoi(argv[4]);

	Tile *tile1 = new Tile(argv[1],n1+1);
	Tile *tile2 = new Tile(argv[2],n2+1);

	int lod = 100;
	if(argc>5){
		lod = atoi(argv[5]);
	}

	tile1->decode_all(lod);
	tile2->decode_all(lod);

	HiMesh *mesh1 = tile1->get_mesh(n1);
	HiMesh *mesh2 = tile2->get_mesh(n2);
	mesh1->get_segments();
	log("%d %d", mesh1->intersect(mesh2), mesh1->intersect_tree(mesh2));

	delete tile1;
	delete tile2;
}

static void print(int argc, char **argv){
	assert(argc>1);
	int num = argc > 2 ? atoi(argv[2]) : 0;
	Tile *tile = new Tile(argv[1],num+1);
	assert(tile->num_objects()>num);
	int lod = argc>3?atoi(argv[3]):100;
	tile->get_mesh(num)->decode(lod);
	cout<<*tile->get_mesh(num);
	delete tile;
}


char *data;
size_t data_size;
int jobs;
int num_jobs;

vector<Tile *> tiles;

pthread_mutex_t mylock;
int next_report = 10;

void *decode_unit(void *arg){
	while(true){
		Tile *tile = NULL;
		int job = -1;
		pthread_mutex_lock(&mylock);
		if(jobs < num_jobs){
			job = jobs++;
			if(jobs*100/num_jobs>=next_report){
				log("%d%%", jobs*100/num_jobs);
				next_report += 10;
			}
		}
		pthread_mutex_unlock(&mylock);
		if(job < 0){
			break;
		}
		tile = tiles[job];

		timeval cur = get_cur_time();
		for(int i=0;i<tile->num_objects();i++){
			HiMesh *mesh = tile->get_mesh(i);
			mesh->decode(100);
		}
		delete tile;
		logt("decode %d", cur, job);
	}
	return NULL;
}

static void decode(int argc, char **argv){

	jobs = 0;
	num_jobs = atoi(argv[2]);

	for(int i=0;i<num_jobs;i++){
		tiles.push_back(new Tile(argv[1]));
	}

	unsigned int num_threads = std::thread::hardware_concurrency();

	if(argc>3){
		num_threads = atoi(argv[3]);
	}

	pthread_t threads[num_threads];

	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, decode_unit, NULL);
	}
	log("%d threads started", num_threads);
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
}

static void print_tile_boxes(int argc, char **argv){
	Tile *tile = new Tile(argv[1]);
	vector<Voxel *> boxes;
	for(int i=0;i<tile->num_objects();i++){
		aab a = tile->get_mesh_wrapper(i)->box;
		Voxel *vx = new Voxel();
		for(int j=0;j<3;j++){
			vx->low[j] = a.low[j];
			vx->high[j] = a.high[j];
		}
		boxes.push_back(vx);
	}
	tdbase::write_voxels(boxes, argv[2]);
}

static void convert(int argc, char **argv){
	Tile *tile = new Tile(argv[1]);
	tile->dump_raw(argv[2]);
	delete tile;
}

static void shrink(int argc, char **argv){
	HiMesh *mesh = read_mesh(argv[1]);
	mesh->shrink(atoi(argv[2]));
	cout<<*mesh;
}

static void shift(int argc, char **argv){
	HiMesh *mesh = read_mesh(argv[1]);
	mesh->shift(mesh->get_mbb().length()/4, mesh->get_mbb().width()/4, mesh->get_mbb().height()/4);
	cout<<*mesh;
}

static void join(int argc, char **argv){
	struct timeval start = get_cur_time();

	geometry_computer *gc = new geometry_computer();
	if(global_ctx.use_gpu){
#ifdef USE_GPU
		initialize();
		gc->init_gpus();
#endif
	}
	if(global_ctx.num_compute_thread>0){
		gc->set_thread_num(global_ctx.num_compute_thread);
	}

	HiMesh::use_byte_coding = !global_ctx.disable_byte_encoding;

	char path1[256];
	char path2[256];

	sprintf(path1, "%s", global_ctx.tile1_path.c_str());
	sprintf(path2, "%s", global_ctx.tile2_path.c_str());

	vector<pair<Tile *, Tile *>> tile_pairs;
	for(int i=0;i<global_ctx.repeated_times;i++){
		Tile *tile1, *tile2;
		if(global_ctx.tile2_path.size()>0){
			tile1 = new Tile(path1, global_ctx.max_num_objects1, false);
			tile2 = new Tile(path2, global_ctx.max_num_objects2, false);
		}else{
			tile1 = new Tile(path1, LONG_MAX, false);
			tile2 = tile1;
		}
		assert(tile1&&tile2);
		tile_pairs.push_back(pair<Tile *, Tile *>(tile1, tile2));
	}
	logt("create tiles", start);

#pragma omp parallel for
	for(int i=0;i<global_ctx.repeated_times;i++){
		Tile *t1 = tile_pairs[i].first;
		Tile *t2 = tile_pairs[i].second;
		t1->load();
		if(t2 != t1){
			t2->load();
		}
	}
	logt("init tiles", start);

	SpatialJoin *joiner = new SpatialJoin(gc);
	joiner->join(tile_pairs);
	double join_time = tdbase::get_time_elapsed(start,false);
	logt("join", start);

#pragma omp parallel for
	for(int i=0;i<global_ctx.repeated_times;i++){
		Tile *t1 = tile_pairs[i].first;
		Tile *t2 = tile_pairs[i].second;
		if(t2 != t1){
			delete t2;
		}
		delete t1;
	}
	tile_pairs.clear();
	logt("clearing tiles", start);
	delete joiner;
	delete gc;
}

static void test(int argc, char **argv){

	tdbase::Point p(0, 1, 2);


}

static void hausdorff(int argc, char **argv){
	HiMesh *low = tdbase::read_mesh(argv[1], false);
	HiMesh *high = tdbase::read_mesh(argv[2], false);
	if(argc>3){
		HiMesh::sampling_rate = atoi(argv[3]);
	}
	low->computeHausdorfDistance(high);
	auto l = low->collectGlobalHausdorff(MIN);
	auto a = low->collectGlobalHausdorff(AVG);
	auto h = low->collectGlobalHausdorff(MAX);

	log("proxy hausdorff [min	avg	max]=[%f	%f	%f]\n",l.first,a.first,h.first);
	log("hausdorff [min	avg	max]=[%f	%f	%f]\n",l.second,a.second,h.second);
}

}

int main(int argc, char **argv){

//	float triangle[] = {248.444000,137.498001,454.556000,
//						252.365005, 133.500000, 456.398987,
//						252.598007, 130.871994, 487.509003};
//
//	float points[][3] = {{248.940002,137.160995,472.493988},
//							{250.354996,136.102005,467.498993},
//							{250.501999,136.369995,477.430511},
//							{251.289001,134.500000,471.625000},
//							{248.432007,132.783005,482.509491},
//							{250.054993,131.158997,487.490997},
//							{250.639008,134.501999,461.420013},
//							{248.501007,138.686005,464.294495},
//							{252.498993,131.270004,456.523010},
//							{252.501007,132.968994,467.226990},
//							{254.479996,125.822998,457.493011},
//							{252.552002,132.505005,480.407990},
//							{254.503006,130.113007,484.035492},
//							{253.501999,129.115997,460.041992},
//							{255.498993,129.488007,473.717499},
//							{253.501007,130.554001,470.134491},
//							{250.089005,134.501999,480.497986},
//							{0,0,0}};
//	for(int i=0;i<18;i++){
//		cout<<points[i][0]<<" "<<points[i][1]<<" "<<points[i][2]<<endl;
//		cout<<tdbase::PointInTriangleCylinder(points[i], triangle)<<endl;
//		//cout<<tdbase::PointTriangleDist(points[i], triangle)<<endl;
//	}
//
//	return 0;

	global_ctx = parse_args(argc, argv);
	if(argc==1){
		cout<<"usage: 3dpro function [args]"<<endl;
		exit(0);
	}

	if(strcmp(argv[1],"join") == 0){
		join(argc-1,argv+1);
	}else if(strcmp(argv[1],"to_wkt") == 0){
		to_wkt(argc-1,argv+1);
	}else if(strcmp(argv[1],"to_sql") == 0){
		to_sql(argc-1,argv+1);
	}else if(strcmp(argv[1],"profile_protruding") == 0){
		profile_protruding(argc-1,argv+1);
	}else if(strcmp(argv[1],"get_voxel_boxes") == 0){
		get_voxel_boxes(argc-1,argv+1);
	}else if(strcmp(argv[1],"profile_distance") == 0){
		profile_distance(argc-1,argv+1);
	}else if(strcmp(argv[1],"profile_decoding") == 0){
		profile_decoding(argc-1,argv+1);
	}else if(strcmp(argv[1],"aabb") == 0){
		aabb(argc-1,argv+1);
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
	}else if(strcmp(argv[1],"simplify") == 0){
		simplify(argc-1,argv+1);
	}else if(strcmp(argv[1],"triangulate") == 0){
		triangulate(argc-1,argv+1);
	}else if(strcmp(argv[1],"distance") == 0){
		distance(argc-1,argv+1);
	}else if(strcmp(argv[1],"print") == 0){
		print(argc-1,argv+1);
	}else if(strcmp(argv[1],"intersect") == 0){
		intersect(argc-1,argv+1);
	}else if(strcmp(argv[1],"print_tile_boxes") == 0){
		print_tile_boxes(argc-1,argv+1);
	}else if(strcmp(argv[1],"sample") == 0){
		sample(argc-1,argv+1);
	}else if(strcmp(argv[1],"decode") == 0){
		decode(argc-1,argv+1);
	}else if(strcmp(argv[1],"convert") == 0){
		convert(argc-1,argv+1);
	}else if(strcmp(argv[1],"shrink") == 0){
		shrink(argc-1,argv+1);
	}else if(strcmp(argv[1],"shift") == 0){
		shift(argc-1,argv+1);
	}else if(strcmp(argv[1],"hausdorff") == 0){
		hausdorff(argc-1,argv+1);
	}else{
		cout<<"usage: 3dpro himesh_to_wkt|profile_protruding|get_voxel_boxes|profile_distance|profile_decoding|adjust_polyhedron|skeleton|voxelize [args]"<<endl;
		exit(0);
	}

	int *abc = new int[1000];
	return 0;
}

