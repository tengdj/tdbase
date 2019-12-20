using namespace std;

/* 
 * RESQUE processing engine v4.0
 *   It supports spatial join and nearest neighbor query with different predicates
 *   1) parseParameters
 *   2) readCacheFile - metadata such as partition schemata
 *   3) for every input line in the current tile
 *         an input line represents an object
 *         save geometry and original data in memory
 *         execute join operation when finish reading a tile
 *   4) Join operation between 2 sets or a single set
 *         build Rtree index on the second data set
 *         for every object in the first data set
 *            using Rtree index of the second data set
 *              check for MBR/envelope intersection
 *              output the pair result or save pair statistics
 *   5) Output final statistics (opt)
 *   Requirement (input files): see the Wiki
 * */

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>

#include <spatialindex/SpatialIndex.h>
#include <boost/algorithm/string/replace.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include "spatial/spatial.h"
#include "iutil.h"
#include "../spatial/himesh.h"
#include "../storage/tile.h"
namespace po = boost::program_options;


/*
 * Struct to hold temporary values during spatial operations
 * for one tile, it contains the mbb, offset, length of objects which
 * are parsed from lines in stdin
 * */
struct query_temp {
	std::string tile_id;
	Tile *tile[2];
};


/* Query operator */
struct query_op {
	Jointype join_predicate = ST_INTERSECTS; /* Join predicate*/
	int decomp_lod = 100; // decompression level/level of details 0 to 100
};

bool extract_params_resque(int argc, char **argv, struct query_op &vars){
	string binpath;
	string jointype;
	try {
		po::options_description desc("Options");
		desc.add_options()
			// resque
			("help,h", "This help message")
			("lod,l", po::value<int>(&vars.decomp_lod) , "Decompression LOD. (0, 100]. Default is 100.")
			("predicate,p", po::value<string>(&jointype)->required(), "Predicate for spatial join and nn queries "
					"[ st_intersects | st_touches | st_crosses | st_contains | st_adjacent | st_disjoint "
					"| st_equals | st_dwithin | st_within | st_overlaps | st_nn_voronoi | st_nn_rtree ] ")
			;
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cerr << desc << endl;
			return false;
		}
		po::notify(vm);


		vars.join_predicate = get_join_predicate(jointype.c_str());
		if(vars.join_predicate==ST_ERROR){
			cerr << desc << endl;
			cerr << "type of predict should be set properly"<<endl;
			exit(-1);
		}
	} catch (exception& e) {
		cerr << "error here: " << e.what() << "\n";
		exit(-1);
	}
	return true;
}


/*
 *
 *
 * for intersection part
 *
 *
 * */




/*
 *
 * CGAL related functions
 *
 * */


typedef MyKernel::Triangle_3                               Triangle;
typedef std::vector<Triangle>                               Triangles;
typedef Triangles::iterator                                   Iterator;

typedef CGAL::Bbox_3                                     Bbox;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;
// for AABB tree distance calculation
typedef CGAL::Triangulation_3<MyKernel> Sc_Triangulation;
typedef Sc_Triangulation::Point        Sc_CGAL_Point;

typedef MyKernel::FT Sc_FT;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Sc_Primitive;
typedef CGAL::AABB_traits<MyKernel, Sc_Primitive> Sc_Traits;
typedef CGAL::AABB_tree<Sc_Traits> Sc_Tree;
typedef Sc_Tree::Point_and_primitive_id Sc_Point_and_primitive_id;

// testing intersecting between polyhedrons
bool intersection_flag = false;
struct Report {

  Report(){}

  // callback functor that reports all truly intersecting triangles
  void operator()(const Box* a, const Box* b) const
  {
    if (intersection_flag) {
    	return;
    }
    if ( ! a->handle()->is_degenerate() && ! b->handle()->is_degenerate()
         && CGAL::do_intersect( *(a->handle()), *(b->handle()))) {
      intersection_flag = true;
    }
  }
};

void get_triangle(HiMesh *P, std::vector<Triangle>& triangles,  std::vector<Box>& boxes, std::vector<Box*>& ptr){
	for ( HiMesh::Facet_const_iterator i = P->facets_begin(); i != P->facets_end(); ++i){
		Triangle t(i->halfedge()->vertex()->point(),
				   i->halfedge()->next()->vertex()->point(),
				   i->halfedge()->next()->next()->vertex()->point());
		triangles.push_back(t);
	}
	// Create the corresponding std::vector of bounding boxes
	for ( Iterator i = triangles.begin(); i != triangles.end(); ++i){
		boxes.push_back(Box(i->bbox(), i));
	}

	for ( std::vector<Box>::iterator i = boxes.begin(); i != boxes.end(); ++i){
		ptr.push_back( &*i);
	}
}

bool intersects(HiMesh *P1, HiMesh *P2) {

	// Use Nef_polyhedron for intersection detection
	Triangles triangles1, triangles2;
	std::vector<Box> boxes1, boxes2;
	std::vector<Box*> boxes1_ptr, boxes2_ptr;

	get_triangle(P1, triangles1, boxes1, boxes1_ptr);
	get_triangle(P2, triangles2, boxes2, boxes2_ptr);

	intersection_flag = false;
	CGAL::box_intersection_d( boxes1_ptr.begin(), boxes1_ptr.end(),
							  boxes2_ptr.begin(), boxes2_ptr.end(),
							  Report());
	return intersection_flag;
}

// performs spatial join on the current tile (bucket)
int join_bucket_intersect(struct query_op &stop, struct query_temp &sttemp) {
	struct timeval start = get_cur_time();

	int pairs = 0; // number of satisfied results
	SpatialIndex::ISpatialIndex *spidx = sttemp.tile[1]->build_rtree();
	logt("building index on set 2", start);

	uint pair_num = 0;
	for (int i = 0; i < sttemp.tile[0]->num_objects(); i++) {
		aab_d b(sttemp.tile[0]->get_mbb(i));
		SpatialIndex::Region r(b.low, b.high, 3);
		MyVisitor vis;
		spidx->intersectsWithQuery(r, vis);
		if(vis.matches.size()==0){
			continue;
		}
		sttemp.tile[0]->decode_to(i,stop.decomp_lod);
		for (uint32_t j = 0; j < vis.matches.size(); j++){
			sttemp.tile[1]->decode_to(vis.matches[j],stop.decomp_lod);
		}
		pair_num += vis.matches.size();
	}
	logt("decoding data for %u pairs", start, pair_num);
	for (int i = 0; i < sttemp.tile[0]->num_objects(); i++) {
		aab_d b(sttemp.tile[0]->get_mbb(i));
		SpatialIndex::Region r(b.low, b.high, 3);
		MyVisitor vis;
		spidx->intersectsWithQuery(r, vis);
		if(vis.matches.size()==0){
			continue;
		}
		for (uint32_t j = 0; j < vis.matches.size(); j++){
			pairs += intersects(sttemp.tile[0]->get_mesh_wrapper(i)->mesh,
					sttemp.tile[1]->get_mesh_wrapper(vis.matches[j])->mesh);
		}
	}
	logt("intersect", start);
	return pairs ;
}


/*
 * in each tile, each large object like vessel is extracted into
 * skeletons. Then a voronoi graph is built on top of those skeletons
 * and objects in data set 1 will get there distance to the nearest object
 * in data set 2 by looking up the voronoi graph
 *
 * */
int join_bucket_nn_voronoi(struct query_op &stop, struct query_temp &sttemp) {

	/* Indicates where original data is mapped to */
	int idx1 = 1;
	int idx2 = 2;
	int nuclei_id = 0;
	double low[3], high[3];  // Temporary value placeholders for MBB
//
//	try {
//		// extract the geometry from dataset2 (compressed blood vessels) and extract skeleton
//		Skeleton *skeleton = NULL;
//
//		std::vector<Delaunay::Point_3> P;
//		vector<struct mbb_3d *> geom_mbb2 = sttemp.mbbdata[idx2];
//		for (int i = 0; i < geom_mbb2.size(); i++) {
//			P.insert(P.begin(), sttemp.skeleton_points[2][i].begin(), sttemp.skeleton_points[2][i].end());
//		}
//		// building their Delaunay triangulation (Voronoi).
//		Delaunay tree(P.begin(), P.end());
//
//		// For each nuclei, find its nearest blood vessel by checking voronoi
//		vector<struct mbb_3d *> geom_mbb1 = sttemp.mbbdata[idx1];
//
//		for (int i = 0; i < geom_mbb1.size(); i++) {
//			struct mbb_3d * env1 = geom_mbb1[i];
//			Delaunay::Point_3 nuclei_centroid((env1->low[0]+env1->high[0])*0.5,
//								  (env1->low[1]+env1->high[1])*0.5,
//								  (env1->low[2]+env1->high[2])*0.5);
//			Delaunay::Point_3 nnp = tree.nearest_vertex(nuclei_centroid)->point();
//			double squared_dist = CGAL::to_double(CGAL::squared_distance(nnp, nuclei_centroid));
//			double distance = sqrt(squared_dist);
//			cout << nuclei_id << TAB << nuclei_centroid.x() << TAB << nuclei_centroid.y()
//					<< TAB << nuclei_centroid.z() << TAB << distance << "\n";
//			nuclei_id++;
//		}
//	} catch (Tools::Exception& e) {
//		std::cerr << "******ERROR******" << std::endl;
//		cerr << e.what() << std::endl;
//		return -1;
//	} // end of catch

	return nuclei_id ;
}

/*
 * perform nearest neighbor query on data sets with AABB tree
 * the kernel for CGAL we used in this module is Simple Cartesian
 * */
int join_bucket_nn_rtree(struct query_op &stop, struct query_temp &sttemp) {

	/* Indicates where original data is mapped to */
	int idx1 = 1;
	int idx2 = 2;
	int pairs = 0; // number of satisfied results
	double low[3], high[3];  // Temporary value placeholders for MBB
	int kneighs = 2; // kNN
	struct timeval start = get_cur_time();
//
//	try {
//		/* Handling for special nearest neighbor query */
//		SpatialIndex::IStorageManager *storage = NULL;
//		SpatialIndex::ISpatialIndex *spidx = NULL;
//		if (! build_index_geoms(sttemp.mbbdata[idx2], spidx, storage)) {
//			return -1;
//		}
//		logt("build rtree", start);
//		unordered_set<int> unique_nn_id2; // the unique set of nearest blood vessels' offset
//		unordered_map<int, vector<int>> nn_id2; // mapping between the unique offset and length
//		vector<Point> nuclei_pts;
//
//		vector<struct mbb_3d *> geom_mbb1 = sttemp.mbbdata[idx1];
//		for (int i = 0; i < geom_mbb1.size(); i++) {
//			struct mbb_3d * env1 = geom_mbb1[i];
//			low[0] = env1->low[0];
//			low[1] = env1->low[1];
//			low[2] = env1->low[2];
//			high[0] = env1->high[0];
//			high[1] = env1->high[1];
//			high[2] = env1->high[2];
//
//			// Temporary value placeholders for MBB
//			double np[3];
//			MyVisitor vis;
//			/* R-tree intersection check */
//
//			np[0] = (low[0]+high[0])*0.5;
//			np[1] = (low[1]+high[1])*0.5;
//			np[2] = (low[2]+high[2])*0.5;
//			SpatialIndex::Point nuclei_centroid(SpatialIndex::Point(np, 3));
//			vis.matches.clear();
//			/* Find kNN objects*/
//			spidx->nearestNeighborQuery(kneighs, nuclei_centroid, vis);
//			for (uint32_t j = 0; j < vis.matches.size(); j++) {
//				// push the offset and length of its nearest blood vessels
//				//long offset = sttemp.offsetdata[idx2][vis.matches[j]], length = sttemp.lengthdata[idx2][vis.matches[j]];
//				nn_id2[i].push_back(vis.matches[j]);
//				// record the unique blood vessels' info
//				unique_nn_id2.insert(vis.matches[j]);
//			}
//			nuclei_pts.push_back(Point(np[0], np[1], np[2]));
//		}
//		if(unique_nn_id2.size()==0){
//			return 0;
//		}
//		logt("checking rtree", start);
//
//		/* for each unique nearest blood vessel, construct the AABB tree*/
//		unordered_map<int, Sc_Tree*> id2_aabbtree; // map between unique id of blood vessel and its AABB tree
//		// for each mentioned object in data set 2, build an AABB tree
//		Sc_Tree *tree = NULL;
//		Polyhedron *geom2[unique_nn_id2.size()];
//		int index = 0;
//		for(auto it = unique_nn_id2.begin(); it != unique_nn_id2.end(); ++it, ++index ){
//			HiMesh *mesh = sttemp.meshes[2][*it];
//			mesh->advance_to(100);
//			geom2[index] = mesh->to_polyhedron();
//		}
//		logt("decoding meshes", start);
//		index = 0;
//		for(auto it = unique_nn_id2.begin(); it != unique_nn_id2.end(); ++it, ++index ){
//			tree = new Sc_Tree(faces(*geom2[index]).first, faces(*geom2[index]).second, *geom2[index]);
//			tree->accelerate_distance_queries();
//			id2_aabbtree[*it] = tree;
//		}
//		logt("build aabb tree", start);
//		/* for each nuclei, calculate distance by searching the AABB tree of its k nearest blood vessels*/
//		for (int j = 0; j < nuclei_pts.size(); j++) {
//			vector<int> ids = nn_id2[j];
//			double min_distance = DBL_MAX;
//			for(int m :ids){
//				Sc_Tree *aabbtree = id2_aabbtree[m];
//				assert(aabbtree!=NULL && "should never happen");
//				Sc_FT sqd = aabbtree->squared_distance(nuclei_pts[j]);
//				double distance = (double)CGAL::to_double(sqd);
//				if(min_distance>distance){
//					min_distance = distance;
//				}
//			}
//			min_distance = sqrt(min_distance);
//			pairs++;
//		}
//		logt("do final computation", start);
//
//		delete spidx;
//		delete storage;
//
//		// release aabb tree of blood vessels
//		for (auto it = id2_aabbtree.begin(); it != id2_aabbtree.end(); ++it ) {
//			delete it->second;
//		}
//		id2_aabbtree.clear();
//		for(Polyhedron *poly:geom2){
//			delete poly;
//			poly = NULL;
//		}
//
//	} catch (Tools::Exception& e) {
//		std::cerr << "******ERROR******" << std::endl;
//		cerr << e.what() << std::endl;
//		return -1;
//	} // end of catch

	return pairs ;
}

/* Release objects in memory (for the current tile/bucket) */
void release_mem(struct query_temp &sttemp) {
	delete sttemp.tile[0];
	delete sttemp.tile[1];
}

/*dispatch joins to target module*/
int join_bucket(struct query_op &stop, struct query_temp &sttemp){
	// Process the current tile in memory
	int pairs = 0;
	switch(stop.join_predicate){
	case ST_NN_VORONOI:
		pairs = join_bucket_nn_voronoi(stop, sttemp);
		break;
	case ST_NN_RTREE:
		pairs = join_bucket_nn_rtree(stop, sttemp);
		break;
	default:
		pairs = join_bucket_intersect(stop, sttemp);
	}
	return pairs;
}


int main(int argc, char** argv)
{
	struct timeval start = get_cur_time();
	struct query_op stop;
	if(!extract_params_resque(argc, argv, stop)){
		return 0;
	}

	Tile *tiles[2];
	tiles[0] = new Tile("nuclei_10k.dt");
	tiles[1] = new Tile("nuclei_10k.dt");
	tiles[0]->retrieve_all();
	tiles[1]->retrieve_all();
	logt("organize data",start);
	query_temp sttemp;
	sttemp.tile[0] = tiles[0];
	sttemp.tile[1] = tiles[1];

	join_bucket(stop, sttemp);
	logt("do join",start);

	delete tiles[0];
	delete tiles[1];
	return 0;
}

/* main body of the engine */
int main1(int argc, char** argv)
{
//	struct query_op stop;
//	if(extract_params_resque(argc, argv, stop)){
//		return 0;;
//	}
//
//	// Processing variables
//	std::string input_line; // Temporary line
//	std::vector<std::string> fields; // Temporary fields
//	int sid = 0; // Join index ID for the current object
//	std::string tile_id = ""; // The current tile_id
//	std::string previd = ""; // the tile_id of the previously read object
//	int tile_counter = 0; // number of processed tiles
//	query_temp sttemp;
//
//	long offset = 0, length = 0;
//	struct mbb_3d *mbb_ptr = NULL;
//
//	// Read line by line inputs
//	while (std::cin && getline(std::cin, input_line) && !std::cin.eof()) {
//
//		tokenize(input_line, fields, TAB, true);
//		if(fields.size()!=11){//skip the corrupted lines
//			std::cerr<<"malformed: "<<fields.size()<<endl;
//						std::cerr<<input_line;
//			continue;
//		}
//		/* Parsing fields from input */
//		tile_id = fields[0];
//		// dataset id
//		sid = atoi(fields[2].c_str());
//		try {
//			// Parsing MBB
//			mbb_ptr = new struct mbb_3d();
//			for (int k = 0; k < 3; k++) {
//				mbb_ptr->low[k] = std::atof(fields[3 + k].c_str());
//			}
//			for (int k = 0; k < 3; k++) {
//				mbb_ptr->high[k] = std::atof(fields[6 + k].c_str());
//			}
//		} catch (...) {
//			std::cerr << "******MBB Parsing Error******" << std::endl;
//			return -1;
//		}
//
//		try {
//			offset = atol(fields[9].c_str());
//			length = atol(fields[10].c_str());
//		}catch (...) {
//			std::cerr << "******Offset and Length Parsing Error******" << std::endl;
//			return -1;
//		}
//
//		/* Process the current tile (bucket) when finishing reading all objects belonging
//		   to the current tile */
//		if (previd.compare(tile_id) != 0 && previd.size() > 0 ) {
//			join_bucket(stop, sttemp);
//			release_mem(sttemp);
//			tile_counter++;
//		}
//
//		// populate the bucket for join
//		sttemp.mbbdata[sid].push_back(mbb_ptr);
//
//		/* Update the field */
//		previd = tile_id;
//		fields.clear();
//	}
//	// Process the last tile (what remains in memory)
//	join_bucket(stop, sttemp);
//	release_mem(sttemp);
//	tile_counter++;

	return 0;
}
