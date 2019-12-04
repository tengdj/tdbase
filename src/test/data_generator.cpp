/*
 * data_generator.cpp
 *
 *  Created on: Nov 27, 2019
 *      Author: teng
 *
 *  generate testing data by duplicating some given examples
 *
 */

#include "../spatial/himesh.h"
#include <fstream>

using namespace hispeed;
using namespace std;


vector<Polyhedron> nucleis;
vector<vector<Voxel *>> nucleis_voxels;
Polyhedron vessel;
vector<Voxel *> vessel_voxels;
aab nuclei_box;
aab vessel_box;
float shrink = 20;
float *possibility;
int num_nuclei = 250;

int num_nuclei_per_vessel(){
	int total_nuclei_num = 1;
	for(int i=0;i<3;i++){
		total_nuclei_num *= (int)(vessel_box.max[i]/nuclei_box.max[i]);
	}
	return total_nuclei_num;
}

HiMesh *poly_to_himesh(Polyhedron &poly){
	stringstream ss;
	ss<<poly;
	MyMesh *mesh = hispeed::get_mesh(ss.str(), true);
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);
	delete mesh;
	return himesh;
}

MyMesh *poly_to_mesh(Polyhedron &poly){
	stringstream ss;
	ss<<poly;
	return hispeed::get_mesh(ss.str(), true);
}

void load_prototype(const char *nuclei_path, const char *vessel_path){
	std::ifstream nfile(nuclei_path);
	std::ifstream vfile(vessel_path);
	string input_line;
	while(std::getline(vfile, input_line)){
		hispeed::replace_bar(input_line);
		stringstream ss;
		ss<<input_line;
		ss >> vessel;
		aab tmpb;
		for(Polyhedron::Vertex_iterator vi=vessel.vertices_begin();vi!=vessel.vertices_end();vi++){
			Point p = vi->point();
			tmpb.update(p[0], p[1], p[2]);
		}
		for(Polyhedron::Vertex_iterator vi=vessel.vertices_begin();vi!=vessel.vertices_end();vi++){
			Point p = vi->point();
			vi->point() = Point(p[0]-tmpb.min[0], p[1]-tmpb.min[1], p[2]-tmpb.min[2]);
		}
		tmpb.max[0] -= tmpb.min[0];
		tmpb.max[1] -= tmpb.min[1];
		tmpb.max[2] -= tmpb.min[2];
		tmpb.min[0] = 0;
		tmpb.min[1] = 0;
		tmpb.min[2] = 0;
		vessel_box.update(tmpb);
		HiMesh *himesh = poly_to_himesh(vessel);
		vessel_voxels = himesh->generate_voxels(500);
		delete himesh;
		break;
	}

	while(std::getline(nfile, input_line)){
		hispeed::replace_bar(input_line);
		stringstream ss;
		ss<<input_line;
		Polyhedron poly;
		ss >> poly;
		aab tmpb;
		for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
			Point p = vi->point();
			Point np(p[0]/shrink, p[1]/shrink, p[2]/shrink);
			vi->point() = np;
			tmpb.update(np[0], np[1], np[2]);
		}
		for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
			Point p = vi->point();
			vi->point() = Point((p[0]-tmpb.min[0]), (p[1]-tmpb.min[1]), p[2]-tmpb.min[2]);
		}
		tmpb.max[0] -= tmpb.min[0];
		tmpb.max[1] -= tmpb.min[1];
		tmpb.max[2] -= tmpb.min[2];
		tmpb.min[2] = 0;
		tmpb.min[1] = 0;
		tmpb.min[2] = 0;
		nuclei_box.update(tmpb);
		nucleis.push_back(poly);
		HiMesh *himesh = poly_to_himesh(poly);
		nucleis_voxels.push_back(himesh->generate_voxels(500));
		delete himesh;
	}

	nfile.close();
	vfile.close();
	// update the distance
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.max[i]/nuclei_box.max[i]);
	}
	int total_nuclei_num = nuclei_num[0]*nuclei_num[1]*nuclei_num[2];
	float pos_adjust = ((float)num_nuclei)/total_nuclei_num*1.5;
	possibility = new float[total_nuclei_num];
	float max_distance = DBL_MIN;
	int zero_num = 0;
	for(int i=0;i<nuclei_num[0];i++){
		for(int j=0;j<nuclei_num[1];j++){
			for(int k=0;k<nuclei_num[2];k++){
				float min_dist = DBL_MAX;
				aab cur_box;
				cur_box.min[0] = nuclei_box.max[0]*i;
				cur_box.max[0] = nuclei_box.max[0]*(i+1);
				cur_box.min[1] = nuclei_box.max[1]*j;
				cur_box.max[1] = nuclei_box.max[1]*(j+1);
				cur_box.min[2] = nuclei_box.max[2]*k;
				cur_box.max[2] = nuclei_box.max[2]*(k+1);
				for(Voxel *v:vessel_voxels){
					range r = v->box.distance(cur_box);
					if(min_dist>r.closest){
						min_dist = r.closest;
						if(min_dist == 0){
							zero_num++;
							break;
						}
					}
				}
				possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k] = min_dist;
				if(max_distance<min_dist){
					max_distance = min_dist;
				}
			}
		}
	}
	// update the possibility for each voxel
	for(int i=0;i<nuclei_num[0];i++){
		for(int j=0;j<nuclei_num[1];j++){
			for(int k=0;k<nuclei_num[2];k++){
				float cur_possibility = possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k];
				if(cur_possibility!=0){
					cur_possibility = (1-cur_possibility/max_distance)*pos_adjust;
					possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k] = cur_possibility;
				}
			}
		}
	}
}


Polyhedron shift_nuclei(float shift[3], int polyid){
	Polyhedron poly;
	stringstream ss;
	ss<<nucleis[polyid];
	ss >> poly;
	for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
		Point p = vi->point();
		vi->point() = Point((p[0]+shift[0]), (p[1]+shift[1]), p[2]+shift[2]);
	}
	return poly;
}

void organize_data(Polyhedron &poly, vector<Voxel *> voxels, float shift[3], char *data, size_t &offset){
	MyMesh *mesh = poly_to_mesh(poly);
	memcpy(data+offset, (char *)&mesh->dataOffset, sizeof(size_t));
	offset += sizeof(size_t);
	memcpy(data+offset, mesh->p_data, mesh->dataOffset);
	offset += mesh->dataOffset;
	size_t size = voxels.size();
	memcpy(data+offset, (char *)&size, sizeof(size_t));
	offset += sizeof(size_t);
	float box_tmp[3];
	for(Voxel *v:voxels){
		for(int i=0;i<3;i++){
			box_tmp[i] = v->box.min[i]+shift[i];
		}
		memcpy(data+offset, (char *)box_tmp, 3*sizeof(float));
		offset += 3*sizeof(float);

		for(int i=0;i<3;i++){
			box_tmp[i] = v->box.max[i]+shift[i];
		}
		memcpy(data+offset, (char *)box_tmp, 3*sizeof(float));
		offset += 3*sizeof(float);

		for(int i=0;i<3;i++){
			box_tmp[i] = v->core[i]+shift[i];
		}
		memcpy(data+offset, (char *)box_tmp, 3*sizeof(float));
		offset += 3*sizeof(float);
	}
}


int generate_nuclei(float base[3], char *data, size_t &offset){
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.max[i]/nuclei_box.max[i]);
	}
	float shift[3];
	int generated = 0;
	for(int i=0;i<nuclei_num[0];i++){
		for(int j=0;j<nuclei_num[1];j++){
			for(int k=0;k<nuclei_num[2];k++){
				float cur_possibility = possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k];
				if(!hispeed::get_rand_sample(cur_possibility)){
					continue;
				}
				shift[0] = i*nuclei_box.max[0]+base[0];
				shift[1] = j*nuclei_box.max[1]+base[1];
				shift[2] = k*nuclei_box.max[2]+base[2];
				int polyid = hispeed::get_rand_number(nucleis.size()-1);
				Polyhedron poly = shift_nuclei(shift, polyid);
				organize_data(poly, nucleis_voxels[polyid], shift, data, offset);
				hispeed::write_polyhedron(&poly, generated);
				generated++;
			}
		}
	}
	return generated;
}


int main(int argc, char **argv){

	load_prototype("../data/nuclei.pt","../data/vessel.pt");
	hispeed::write_polyhedron(&vessel, "v.off");
	char *data = new char[1<<30];
	size_t offset = 0;
	float base[3] = {0,0,0};
	int generated = generate_nuclei(base, data, offset);
	cout<<generated<<endl;
	ofstream *os = new std::ofstream("nuclei_tmp.dt", std::ios::out | std::ios::binary);
	os->write(data, offset);
	os->close();
	delete os;


}

