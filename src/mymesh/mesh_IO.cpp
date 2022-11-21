/*
 * mesh_parse.cpp
 *
 *  Created on: Nov 16, 2022
 *      Author: teng
 */

#include "mymesh.h"

namespace hispeed{


string Polyhedron::to_string(){

	stringstream os;
	os << "OFF" << endl;
	os << vertices.size() << " " << faces.size() << " 0\n" << endl;
	int idx = 0;
	for(Vertex *p:vertices){
		p->id = idx++;
		os << p->v[0] << " " << p->v[1] << " " << p->v[2] << endl;
	}
	for(Face *f:faces){
		os<<f->vertices.size()<<" ";
		for(Vertex *v:f->vertices){
			os<<v->id<<" ";
		}
		if(f->added){
			os<<"\t255 255 0";
		}else{
			os<<"\t0 255 0";
		}
		os<<endl;
	}
	return os.str();
}

void Polyhedron::print(){
	cout<<to_string().c_str();
}

void Polyhedron::dumpto(string fp){
	ofstream of(fp.c_str());
	of<<to_string().c_str();
	of.close();
}

//
// read vertices and faces from original OFF file
//
void Polyhedron::load(string fp) {
	string content = hispeed::read_file(fp.c_str());
	parse(content.c_str(), content.size());
}


bool parse_OFF_sign(const char *data, size_t &offset, size_t len){
	if(len<3){
		return false;
	}
	while(offset<len-3){
		if((data[offset]=='O' || data[offset]=='o')
			&&(data[offset+1]=='F' || data[offset+1]=='f')
			&&(data[offset+2]=='F' || data[offset+2]=='f')){
			offset += 3;
			return true;
		}else{
			offset++;
			continue;
		}
	}
	// no OFF is found;
	return false;
}

void skip_spaces(const char *data, size_t &offset, size_t len){
	while(offset<len){
		if(data[offset]==' '||data[offset]=='\t'||data[offset]=='\n'||data[offset]=='|'){
			offset++;
			continue;
		}
		break;
	}
}

bool parse_float(const char *data, size_t &offset, size_t len, float &ret){
	skip_spaces(data, offset, len);
	char buffer[20];
	int iter = 0;
	bool has_point = false;
	while(true){
		if(offset==len || ((data[offset]<'0' || data[offset]>'9') && data[offset]!='.')){
			if(iter>0){
				buffer[iter] = '\0';
				ret = atof(buffer);
				return true;
			}else{
				return false;
			}
		}else{
			buffer[iter++] = data[offset++];
			if(iter>=19){
				buffer[iter] = '\0';
				log("%s is too long", buffer);
				return false;
			}

			if(data[offset]=='.'){
				if(has_point){
					buffer[iter] = '\0';
					log("wrong float value", buffer);
					return false;
				}
				has_point = true;
			}
		}
	}
	return false;
}

bool parse_int(const char *data, size_t &offset, size_t len, int &ret){
	skip_spaces(data, offset, len);
	char buffer[20];
	int iter = 0;
	while(true){
		if(offset==len || data[offset]<'0' || data[offset]>'9'){
			if(iter>0){
				buffer[iter] = '\0';
				ret = atoi(buffer);
				return true;
			}else{
				return false;
			}
		}else{
			buffer[iter++] = data[offset++];
			if(iter>=19){
				buffer[iter] = '\0';
				log("%s is too long", buffer);
				return false;
			}
		}
	}
	return false;
}

bool Polyhedron::parse(string str){
	return parse(str.c_str(), str.size());
}

// parse from OFF
bool Polyhedron::parse(const char *data, size_t len){
	size_t offset = 0;
	float flt_tmp;
	int int_tmp;
	if(!parse_OFF_sign(data,offset,len)){
		log("no OFF sign found");
		return false;
	}
	int num_vertices, num_faces, num_edges;
	if(!parse_int(data, offset, len, num_vertices)||
	   !parse_int(data, offset, len, num_faces)||
	   !parse_int(data, offset, len, num_edges)){
		return false;
	}
	if(num_vertices==0){
		return true;
	}
	vector<Vertex *> vertex_array;
	vertex_array.resize(num_vertices);
	for(int i=0;i<num_vertices;i++){
		float co[3];
		for(int j=0;j<3;j++){
			if(!parse_float(data, offset, len, co[j])){
				return false;
			}
		}
		vertex_array[i] = new Vertex(co[0], co[1], co[2]);
	}
	if(num_faces==0){
		return true;
	}
	for(int i=0;i<num_faces;i++){
		int fnum = 0;
		if(!parse_int(data, offset, len, fnum) || fnum<3){
			log("invalid face numbers %d", fnum);
			return false;
		}
		vector<Vertex *> vs;
		vs.resize(fnum);
		for(int j=0;j<fnum;j++){
			int vindex = 0;
			if(!parse_int(data, offset, len, vindex) || vindex>=vertex_array.size()){
				log("invalid vertex index %d", vindex);
				return false;
			}
			vs[j] = vertex_array[vindex];
		}
		add_face(vs);
	}
	vertices.insert(vertex_array.begin(), vertex_array.end());

	return true;
}


}


