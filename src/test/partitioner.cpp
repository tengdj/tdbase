/*
 * partitioner.cpp
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 */

#include "../partition/partition.h"

using namespace hispeed;
using namespace std;

// main method
int main(int argc, char** argv) {


	std::vector<string> input_folders;
	//input_folders.push_back("/home/teng/project/HiSPEED/build/nuclei");
	//input_folders.push_back("/home/teng/project/HiSPEED/build/vessels");
	input_folders.push_back(argv[1]);
	std::vector<aab> tiles;
	hispeed::partition_space(input_folders, tiles, 10, 100);
	long sum = 0;
	for(int i=0;i<tiles.size();i++){
		cout<<i<<" "<<tiles[i]<<endl;
		sum += tiles[i].weight;
	}
	cout<<tile_size*100<<" "<<sum<<endl;
	return 0;
}



