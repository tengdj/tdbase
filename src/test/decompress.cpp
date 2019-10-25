/*
 * decomp.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: teng
 */


#include "../PPMC/ppmc.h"
#include "../util/util.h"

int main(int argc, char **argv){

	int start_lod = 0;
	int end_lod = 10;
	if(argc>1){
		start_lod = atoi(argv[1]);
		end_lod = start_lod;
	}
	assert(start_lod>=0&&start_lod<=10);

	int i_mode = COMPRESSION_MODE_ID; // compression mode
	unsigned i_quantBit = 12;
	unsigned i_decompPercentage = 100;
	bool b_allowConcaveFaces = true;
	string input_line;
	getline(std::cin, input_line);
	boost::replace_all(input_line, "|", "\n");
	// Init the random number generator.
	std::cerr<<"start compressing"<<endl;
	srand(PPMC_RANDOM_CONSTANT);
	struct timeval starttime = get_cur_time();
	MyMesh *compressed = new MyMesh(i_decompPercentage,
				 COMPRESSION_MODE_ID, i_quantBit,
				 b_allowConcaveFaces,
				 input_line.c_str(), input_line.size());
	compressed->completeOperation();
	cout<<compressed->dataOffset<<" "<<get_time_elapsed(starttime)<<endl;

	std::cerr<<"start decompressing"<<endl;
	char path[100];

	for(int i=start_lod;i<=end_lod;i++){
		starttime = get_cur_time();
		srand(PPMC_RANDOM_CONSTANT);
		MyMesh *decompressed = new MyMesh(10*i,
				 DECOMPRESSION_MODE_ID, i_quantBit,
				 b_allowConcaveFaces,
				 compressed->p_data, compressed->dataOffset);
		decompressed->completeOperation();
		sprintf(path,"lod%d.off", i*10);
		print_mesh_file(decompressed, path);
		decompressed->teng_test();

		delete decompressed;
		cout<<decompressed->dataOffset<<" "<<get_time_elapsed(starttime)<<endl;

	}

	delete compressed;

}



