/*
 * 3dpro_tool.cpp
 *
 *  Created on: Mar 29, 2022
 *      Author: teng
 */

#include <algorithm>
#include <zlib.h>
#include <thread>

#include "himesh.h"
#include "tile.h"
#include "util.h"

using namespace std;
using namespace hispeed;

namespace hispeed{
/*
 * print himesh to wkt
 *
 * */
static void himesh_to_wkt(int argc, char **argv){

	Tile *tile = new Tile(argv[1]);
	const char *prefix = argv[2];
	tile->retrieve_all();
	tile->advance_all(100);
	char path[256];
	for(int i=0;i<tile->num_objects();i++){
		sprintf(path, "%s_%d.OFF", prefix, i);
		tile->get_mesh(i)->write_to_off(path);
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
	HiMesh *mesh = hispeed::read_mesh();
	int voxel_num = 100;
	if(argc>=2){
		voxel_num = atoi(argv[1]);
	}
	vector<Voxel *> voxels = mesh->generate_voxels_skeleton(voxel_num);

	float vol = 0.0;
	for(Voxel *v:voxels){
		vol += v->volume();
	}

	log("%ld voxels are generated %f volumn", voxels.size(), vol);
	write_voxels(voxels, "/gisdata/skeleton_voxels.off");
	hispeed::write_box(mesh->get_mbb(), "/gisdata/aab.off");
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
	vector<Point> skeleton = mesh->get_skeleton_points(voxel_num);
	hispeed::write_points(skeleton, argv[2]);

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
		tile1->advance_all(i);
		tile2->advance_all(i);
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
	hispeed::write_box(box, "/gisdata/box.off");

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

	global_ctx.use_multimbb = false;
	Tile *tile1 = new Tile(argv[1], num);
	Tile *tile2 = new Tile(argv[2],1);
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
	HiMesh *mesh = hispeed::read_mesh();
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
	hispeed::write_polyhedron(poly, argv[2]);

	delete poly;
	delete mesh;
}

static void sample(int argc, char **argv){

	HiMesh::sampling_rate = atoi(argv[2]);
	int lod = 100;
	if(argc>3){
		lod = atoi(argv[3]);
	}

	log("sampling rate %d lod %d", HiMesh::sampling_rate, lod);

	struct timeval start = get_cur_time();
	HiMesh *final_mesh;
	if(lod!=100){
		HiMesh *mesh = read_mesh(argv[1], true);
		logt("compress", start);
		final_mesh = new HiMesh(mesh);
		final_mesh->decode(lod);
		logt("decompress", start);
		delete mesh;
	}else{
		final_mesh = read_mesh(argv[1], false);
	}

	char path[256];

	unordered_set<Point> point_set;
	final_mesh->sample_points(point_set);
	vector<Point> points;
	points.assign(point_set.begin(), point_set.end());
    sprintf(path, "/gisdata/sample.points.off");
	hispeed::write_points(points, path);
    sprintf(path, "/gisdata/sample.mesh.off");
    final_mesh->write_to_off(path);

	delete final_mesh;
}

static void compress(int argc, char **argv){

//	Point p1(201.502, 160.057, 245.023);
//	Point p2(221.539, 150.502, 155.167);
//	Point p3(206.559, 152.501, 73.854);
//
//	cout<<triangle_area(p1, p2, p3)<<endl;
//
//	return;

	global_ctx.verbose = 2;
	if(argc>2){
		HiMesh::sampling_rate = atoi(argv[2]);
		log("%d",HiMesh::sampling_rate);
	}

	if(argc>3){
		HiMesh::calculate_method = atoi(argv[3]);
	}
	struct timeval start = get_cur_time();
	HiMesh *mesh = read_mesh(argv[1], true);
	logt("compress", start);
	HiMesh *hm = new HiMesh(mesh);
	int lod = 100;

	char path[256];
	for(uint i=0;i<=lod;i+=10){
		hm->decode(i);
		logt("decode to %d", start, i);
		//log("%d %f", i, himesh->getHausdorfDistance());
	    sprintf(path, "/gisdata/compressed_%d.mesh.off", i);
	    hm->write_to_off(path);

//		unordered_set<Point> point_set;
//		hm->sample_points(point_set);
//		vector<Point> points;
//		points.assign(point_set.begin(), point_set.end());
//	    sprintf(path, "/gisdata/compressed_%d.points.off", i);
//		hispeed::write_points(points, path);
//
//		//log("global: %f", himesh->getHausdorfDistance().second);
//		int tris = hm->size_of_triangles();
//		for(int j=0;j<tris;j++){
//			//log("%d\t%.2f\t%d", j, himesh->get_triangle_hausdorf(j).second, (int)(himesh->get_triangle_hausdorf(j).second*100/himesh->getHausdorfDistance().second));
//		}
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

	tile1->advance_all(lod);
	tile2->advance_all(lod);

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

	tile1->advance_all(lod);
	tile2->advance_all(lod);

	HiMesh *mesh1 = tile1->get_mesh(n1);
	HiMesh *mesh2 = tile2->get_mesh(n2);
	mesh1->get_segments();
	log("%d %d", mesh1->intersect(mesh2), mesh1->intersect_tree(mesh2));

	delete tile1;
	delete tile2;
}

static void print(int argc, char **argv){
	assert(argc>2);
	Tile *tile = new Tile(argv[1],atoi(argv[2])+1);
	assert(tile->num_objects()>atoi(argv[2]));
	int lod = argc>3?atoi(argv[3]):100;
	tile->get_mesh(atoi(argv[2]))->decode(lod);
	cout<<*tile->get_mesh(atoi(argv[2]));
	delete tile;
}


char *data;
size_t data_size;
int jobs;
int num_jobs;

vector<Tile *> tiles;

pthread_mutex_t mylock;
int next_report = 10;

void *generate_unit(void *arg){
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

#ifdef CGAL_HAS_THREADS
	log("terry is good");
#endif

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
		pthread_create(&threads[i], NULL, generate_unit, NULL);
	}
	log("%d threads started", num_threads);
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
}

//void *generate_unit(void *arg){
//	while(true){
//		int job;
//		Tile *tile = NULL;
//		pthread_mutex_lock(&mylock);
//		job = jobs--;
//		if(job>=0 && (num_jobs-job)*100/num_jobs>=next_report){
//			log("%d%%", (num_jobs-job)*100/num_jobs);
//			next_report += 10;
//		}
//		pthread_mutex_unlock(&mylock);
//		if(job<0){
//			break;
//		}
//		timeval cur = get_cur_time();
//		HiMesh *mesh = new HiMesh(data, data_size);
//		mesh->decode(100);
//		delete mesh;
//		//logt("decode %d", cur, job);
//
//	}
//	return NULL;
//}
//
//static void decode(int argc, char **argv){
//	Tile *tile = new Tile(argv[1], 1);
//	data = tile->retrieve_data(0);
//	data_size = tile->get_object_data_size(0);
//
//	num_jobs = atoi(argv[2]);
//	jobs = num_jobs;
//
//	unsigned int num_threads = std::thread::hardware_concurrency();
//
//	if(argc>3){
//		num_threads = atoi(argv[3]);
//	}
//
//	pthread_t threads[num_threads];
//
//	for(int i=0;i<num_threads;i++){
//		pthread_create(&threads[i], NULL, generate_unit, NULL);
//	}
//	log("%d threads started", num_threads);
//	for(int i = 0; i < num_threads; i++ ){
//		void *status;
//		pthread_join(threads[i], &status);
//	}
//}

//static void decode(int argc, char **argv){
//	Tile *tile = new Tile(argv[1], 1);
//	char *data = tile->retrieve_data(0);
//	size_t data_size = tile->get_object_data_size(0);
//
//	size_t no = atoi(argv[2]);
//
//#pragma omp parallel for
//	for(size_t i=0;i<no;i++){
//		struct timeval start = get_cur_time();
//		HiMesh *mesh = new HiMesh(data, data_size);
//		mesh->decode(100);
//		logt("decoding %ld", start, i);
//	}
//
//}

//static void decode(int argc, char **argv){
//
//	int nt = atoi(argv[2]);
//#pragma omp parallel for
//	for(int t=0;t<nt;t++){
//		Tile *tile = new Tile(argv[1]);
//		for(size_t i=0;i<tile->num_objects();i++){
//			struct timeval start = get_cur_time();
//			HiMesh *mesh = tile->get_mesh(i);
//			mesh->decode(100);
//			logt("decoding %d-%ld", start, t, i);
//		}
//		delete tile;
//	}
//
//}

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
	hispeed::write_voxels(boxes, argv[2]);
}

//static replacing_group *merge1(unordered_set<replacing_group *> &reps){
//	replacing_group *ret = new replacing_group();
//	for(replacing_group *r:reps){
//		ret->removed_vertices.insert(r->removed_vertices.begin(), r->removed_vertices.end());
//		delete r;
//		r = NULL;
//	}
//	for(Point p:ret->removed_vertices){
//		cout<<p<<endl;
//	}
//	log("merged %ld reps with %ld removed vertices", reps.size(), ret->removed_vertices.size());
//	reps.clear();
//	return ret;
//}

static void test(int argc, char **argv){
//
//	Point p1(1,2,3);
//	Point p2(2,3,4);
//	Point p3(3,4,5);
//	Point p4(4,5,6);
//
//	replacing_group *r1 = new replacing_group();
//	r1->removed_vertices.emplace(p1);
//	r1->removed_vertices.emplace(p2);
//
//	replacing_group *r2 = new replacing_group();
//	r2->removed_vertices.emplace(p3);
//	r2->removed_vertices.emplace(p4);
//
//	unordered_set<replacing_group *> reps;
//	reps.emplace(r1);
//	reps.emplace(r2);
//	auto r3 = merge1(reps);
//	for(auto a:r3->removed_vertices){
//		cout<<a<<endl;
//	}

//	Tile *tile = new Tile(argv[1],1);
//	HiMesh *mesh = tile->get_mesh(0);
//
//	cout<<sizeof(Polyhedron)<<endl;
//
//	if(argc>2){
//		int lod = atoi(argv[2]);
//		assert(lod>=0 && lod<=100);
//		mesh->decode(lod);
//	}
//	cout<<*mesh<<endl;

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
//		cout<<hispeed::PointInTriangleCylinder(points[i], triangle)<<endl;
//		//cout<<hispeed::PointTriangleDist(points[i], triangle)<<endl;
//	}
//
//	return 0;

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
	}else{
		cout<<"usage: 3dpro himesh_to_wkt|profile_protruding|get_voxel_boxes|profile_distance|profile_decoding|adjust_polyhedron|skeleton|voxelize [args]"<<endl;
		exit(0);
	}
	return 0;
}

