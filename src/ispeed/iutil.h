#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <string>
#include <vector>
#include "../util/util.h"
using namespace std;
using namespace hispeed;

// name of binary tools
const std::string MANIPULATE = "manipulate_3d";
const std::string RESQUE = "resque_3d";

#define SPACE " "
#define SLASH "/"

enum Jointype{
	ST_ERROR = 0,
	ST_INTERSECTS = 1,
	ST_TOUCHES = 2,
	ST_CROSSES = 3,
	ST_CONTAINS = 4,
	ST_ADJACENT = 5,
	ST_DISJOINT = 6,
	ST_EQUALS = 7,
	ST_DWITHIN = 8,
	ST_WITHIN = 9,
	ST_OVERLAPS = 10,
	ST_NEAREST = 11,
	ST_NEAREST_2 = 12,
	ST_NN_VORONOI = 13,
	ST_NN_RTREE = 14
};

const std::string join_type_str[15] = {
		"st_error", "st_intersects", "st_touches", "st_crosses", "st_contains",
		"st_adjacent", "st_disjoint", "st_equals", "st_dwithin", "st_within",
		"st_overlaps", "st_nearest", "st_nearest2", "st_nn_voronoi", "st_nn_rtree"};

inline Jointype get_join_predicate(const char * predicate_str){
	for(int i=1;i<15;i++){
		if (strcmp(predicate_str, join_type_str[i].c_str()) == 0) {
			return (Jointype)i ;
		}
	}
	std::cerr << "unrecognized join predicate " <<predicate_str<< std::endl;
	return ST_ERROR;
}

#endif

