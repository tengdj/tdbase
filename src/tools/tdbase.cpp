/*
 * 3dpro_tool.cpp
 *
 *  Created on: Mar 29, 2022
 *      Author: teng
 */

#include <algorithm>
#include <thread>
#include <map>

#include "SpatialJoin.h"
#include "himesh.h"
#include "tile.h"
#include "util.h"

using namespace std;
using namespace tdbase;

namespace tdbase{

	map<string, void (*)(int, char **)> functions;

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

	vector<Point> points;
	points.assign(mesh->sampled_points.begin(), mesh->sampled_points.end());
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

/*for evaluating the decoding efficiency*/
static void decode(int argc, char **argv){
	vector<Tile*> tiles;
	struct timeval start = get_cur_time();
	Tile* tile = new Tile(argv[1]);
	logt("loading",start);
	int num_thread = argc > 2? atoi(argv[2]):tdbase::get_num_threads();
#pragma omp parallel for num_threads(num_thread)
	for (int i = 0; i < tile->num_objects(); i++) {
		HiMesh* mesh = tile->get_mesh(i);
		mesh->decode(100);
	}
	delete tile;
	logt("decoding", start);
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

static void hausdorff(int argc, char** argv) {
	HiMesh* low = tdbase::read_mesh(argv[1], false);
	HiMesh* high = tdbase::read_mesh(argv[2], false);
	if (argc > 3) {
		HiMesh::sampling_rate = atoi(argv[3]);
	}
	low->computeHausdorfDistance(high);
	auto l = low->collectGlobalHausdorff(MIN);
	auto a = low->collectGlobalHausdorff(AVG);
	auto h = low->collectGlobalHausdorff(MAX);

	log("proxy hausdorff [min	avg	max]=[%f	%f	%f]\n", l.first, a.first, h.first);
	log("hausdorff [min	avg	max]=[%f	%f	%f]\n", l.second, a.second, h.second);
}

static void join(int argc, char **argv){
	global_ctx = parse_args(argc, argv);

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

// 自定义 pair<int,int> 的 hash
struct PairHash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<long long>()(
            (static_cast<long long>(p.first) << 32) | (unsigned int)p.second
        );
    }
};

unordered_set<pair<int, int>, PairHash>
load_pairs(const string& filename) {
    unordered_set<pair<int, int>, PairHash> s;
    ifstream fin(filename);
    if (!fin) {
        cerr << "Failed to open file: " << filename << endl;
        exit(1);
    }

    int a, b;
    while (fin >> a >> b) {
        s.emplace(a, b);
    }
    return s;
}

static void evaluate(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " ground_truth.txt prediction.txt\n";
        return;
    }

    auto gt   = load_pairs(argv[1]);
    auto pred = load_pairs(argv[2]);

    size_t true_positive = 0;
    for (const auto& p : pred) {
        if (gt.count(p)) {
            true_positive++;
        }
    }

    double precision = pred.empty() ? 0.0 :
        static_cast<double>(true_positive) / pred.size();
    double recall = gt.empty() ? 0.0 :
        static_cast<double>(true_positive) / gt.size();

    cout << "True Positive: " << true_positive << endl;
    cout << "Precision: " << precision << endl;
    cout << "Recall: " << recall << endl;
}

}

int main(int argc, char **argv){
	// register the functions
	functions["to_sql"] = to_sql;
	functions["get_voxel_boxes"] = get_voxel_boxes;
	functions["get_skeleton"] = get_skeleton;
	functions["voxelize"] = voxelize;
	functions["profile_decoding"] = profile_decoding;
	functions["aabb"] = aabb;
	functions["adjust_polyhedron"] = adjust_polyhedron;
	functions["triangulate"] = triangulate;
	functions["sample"] = sample;
	functions["simplify"] = simplify;
	functions["compress"] = compress;
	functions["distance"] = distance;
	functions["intersect"] = intersect;
	functions["print"] = print;
	functions["decode"] = decode;
	functions["print_tile_boxes"] = print_tile_boxes;
	functions["convert"] = convert;
	functions["shrink"] = shrink;
	functions["shift"] = shift;
	functions["hausdorff"] = hausdorff;
	functions["join"] = join;
	functions["evaluate"] = evaluate;

	if (argc < 2 || functions.find(argv[1])==functions.end()) {
		cout<<"usage: tdbase function [args]"<<endl;
		cout << "function could be:"<<endl;
		for (auto e:functions) {
			cout << e.first << " ";
		}
		cout << endl;
		exit(0);
	}
	functions[argv[1]](argc-1, argv+1);

	return 0;
}

