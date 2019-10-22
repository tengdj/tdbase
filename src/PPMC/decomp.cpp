/*
 * decomp.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: teng
 */


#include "ppmc.h"

using namespace std;



void print_mesh(MyMesh *mesh){
	std::stringstream os;
	os << *mesh;
	cout << os.str()<<endl;
}

void print_mesh_file(MyMesh *mesh, char *path){
	ofstream myfile;
	myfile.open(path);
	myfile << *mesh;
	myfile.close();
}


int main(int argc, char **argv){

	int i_mode = COMPRESSION_MODE_ID; // compression mode
	unsigned i_quantBit = 12;
	unsigned i_decompPercentage = 100;
	bool optimization = false;
	bool b_allowConcaveFaces = true;
	string input_line;
	getline(std::cin, input_line);
	boost::replace_all(input_line, "|", "\n");
	// Init the random number generator.
	std::cerr<<"start compressing"<<endl;
	srand(PPMC_RANDOM_CONSTANT);
	MyMesh *compressed = new MyMesh(i_decompPercentage,
				 COMPRESSION_MODE_ID, i_quantBit,
				 b_allowConcaveFaces,
				 input_line.c_str(), input_line.size());
	compressed->completeOperation();

	std::cerr<<"start decompressing"<<endl;
	char path[100];
	for(int i=0;i<=10;i++){
		srand(PPMC_RANDOM_CONSTANT);
		MyMesh *decompressed = new MyMesh(10*i,
				 DECOMPRESSION_MODE_ID, i_quantBit,
				 b_allowConcaveFaces,
				 compressed->p_data, compressed->dataOffset);
		decompressed->completeOperation();
		sprintf(path,"lod%d.off", i*10);
		print_mesh_file(decompressed, path);
		cout<<decompressed->dataOffset<<endl;
		delete decompressed;
	}

	delete compressed;

}



