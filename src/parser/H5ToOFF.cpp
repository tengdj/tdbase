#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"


class MyPolyhedron{
public:
	MyPolyhedron(){};
	MyPolyhedron(const char *path);
	~MyPolyhedron(){
		if(vertices!=NULL){
			free(vertices[0]);
			free(vertices);
		}
		if(faces!=NULL){
			free(faces[0]);
			free(faces);
		}
	}
	hsize_t vertex_num = 0;
	hsize_t face_num = 0;
	double **vertices = NULL;
	long **faces = NULL;
	void dump_to_off(const char *path);
};

void MyPolyhedron::dump_to_off(const char *path){

	FILE * pFile = fopen(path,"w");

	fprintf(pFile, "OFF\n%lld %lld 0\n\n",vertex_num,face_num);
	for(hsize_t i=0;i<vertex_num;i++){
		fprintf(pFile, "%.1f %.1f %.1f\n",vertices[i][0],vertices[i][1],vertices[i][2]);
	}
	for(hsize_t i=0;i<face_num;i++){
		fprintf(pFile, "3	%ld %ld %ld\n",faces[i][0],faces[i][1],faces[i][2]);
	}

	fclose(pFile);
}

MyPolyhedron::MyPolyhedron(const char *path){
	hid_t	file, space, dset;
	herr_t	status;
	int		ndims;
	hsize_t dims[2];

	file = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT);

	/*
	* Get dataspace and allocate memory for read buffer.  This is a
	* two dimensional dataset so the dynamic allocation must be done
	* in steps.
	*/

	// load vertices
	dset = H5Dopen (file, "vertices", H5P_DEFAULT);
	space = H5Dget_space (dset);
	ndims = H5Sget_simple_extent_dims (space, dims, NULL);

	vertices = (double **) malloc (dims[0] * sizeof (double *));
	vertices[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
	for (int i=1; i<dims[0]; i++)
		vertices[i] = vertices[0] + i * dims[1];
	vertex_num = dims[0];
	//read data
	status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vertices[0]);
	status = H5Dclose (dset);
	status = H5Sclose (space);

	// load faces
	dset = H5Dopen (file, "faces", H5P_DEFAULT);
	space = H5Dget_space (dset);
	ndims = H5Sget_simple_extent_dims (space, dims, NULL);

	faces = (long **) malloc (dims[0] * sizeof (long *));
	faces[0] = (long *) malloc (dims[0] * dims[1] * sizeof (long));
	for (int i=1; i<dims[0]; i++)
		faces[i] = faces[0] + i * dims[1];
	face_num = dims[0];

	//read data
	status = H5Dread (dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces[0]);
	status = H5Dclose (dset);
	status = H5Sclose (space);

	status = H5Fclose (file);
}

void process_file(const char *path){
    int len = strlen(path);

    if(len<=3||path[len-1]!='5'
    		||path[len-2]!='h'||path[len-3]!='.'){
    	fprintf(stderr, "%s is not a HDF5 file\n", path);
    	return;
    }

    char outpath[256];
    strcpy(outpath, path);
    outpath[len-2]='O';
    outpath[len-1]='F';
    outpath[len]='F';
    outpath[len+1]='\0';
    printf("processing %s\n",path);

	// load Polyhedron
    MyPolyhedron *poly = new MyPolyhedron(path);
	poly->dump_to_off(outpath);
	delete poly;
}

//int main(int argc, char **argv){
//
//	vector<string> files;
//	hispeed::list_files(argv[1], files);
//	for(string s:files){
//		process_file(s.c_str());
//	}
//	return 0;
//}
