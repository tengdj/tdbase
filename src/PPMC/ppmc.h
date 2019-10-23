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

inline struct timeval get_cur_time(){
	struct timeval t1;
	gettimeofday(&t1, NULL);
	return t1;
}
inline double get_time_elapsed(struct timeval t1){
	struct timeval t2;
    double elapsedTime;
	gettimeofday(&t2, NULL);
	// compute and print the elapsed time in millisec
	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
	return elapsedTime;
}
#endif
