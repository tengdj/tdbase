#ifndef PPMC_H_
#define PPMC_H_
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <stdio.h>
#include <assert.h>
#include <string.h>


#include <vector>
#include <istream>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <GL/glut.h>
#include <CGAL/Timer.h>

// Program parameters
#include "configuration.h"
#include "mymesh.h"

#define BAR "|"
#define TAB "\t"
#define PPMC_RANDOM_CONSTANT 315


using namespace std;
using namespace CGAL;

inline void print_mesh(MyMesh *mesh){
	std::stringstream os;
	os << *mesh;
	cout << os.str()<<endl;
}

inline void print_mesh_file(MyMesh *mesh, char *path){
	ofstream myfile;
	myfile.open(path);
	myfile << *mesh;
	myfile.close();
}

inline MyMesh *read_mesh(){
	unsigned i_decompPercentage = 100;
	bool b_allowConcaveFaces = true;
	string input_line;
	getline(std::cin, input_line);
	if(input_line.size()==0){
		return NULL;
	}
	boost::replace_all(input_line, "|", "\n");
	// Init the random number generator.
	srand(PPMC_RANDOM_CONSTANT);
	MyMesh *compressed = new MyMesh(i_decompPercentage,
				 COMPRESSION_MODE_ID, 12,
				 b_allowConcaveFaces,
				 input_line.c_str(), input_line.size());
	compressed->completeOperation();
	return compressed;
}


inline MyMesh *decompress_mesh(MyMesh *compressed, int lod){
	MyMesh *decompressed = new MyMesh(lod,
			 DECOMPRESSION_MODE_ID, 12, true,
			 compressed->p_data, compressed->dataOffset);
	decompressed->completeOperation();
	return decompressed;
}


#endif
