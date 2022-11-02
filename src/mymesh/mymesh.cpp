/*
 * mymesh.cpp
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */

#include "mymesh.h"

namespace hispeed{


Polyhedron::~Polyhedron(){
	for(Point *p:points){
		delete p;
	}
	for(Face *f:faces){
		delete f;
	}
	points.clear();
	faces.clear();
}


void Polyhedron::print(){
	cout << "OFF" << endl;
	cout << points.size() << " " << faces.size() << " 0\n" << endl;
	for(Point *p:points){
		cout << p->v[0] << " " << p->v[1] << " " << p->v[2] << endl;
	}
	for(Face *f:faces){
		cout<<f->vertices.size()<<"\t";
		for(int v:f->vertices){
			cout<<v<<" ";
		}
		cout<<endl;
	}
}

void Polyhedron::dumpto(string fp){

	ofstream of(fp.c_str());
	of << "OFF" << endl;
	of << points.size() << " " << faces.size() << " 0\n" << endl;
	for(Point *p:points){
		of << p->v[0] << " " << p->v[1] << " " << p->v[2] << endl;
	}
	for(Face *f:faces){
		of<<f->vertices.size()<<"\t";
		for(int v:f->vertices){
			of<<v<<" ";
		}
		of<<endl;
	}
	of.close();
}

//
// read vertices and faces from original OFF file
//
void Polyhedron::load(string fp) {
    ifstream infile(fp);
    int n_p, n_face, int_temp;
    string line;

    getline(infile, line);  // OFF
    infile >> n_p >> n_face >> int_temp;
    getline(infile, line); // empty_line
    double dv1, dv2, dv3;
    for (int i = 0; i < n_p; ++i) {
        infile >> dv1 >> dv2 >> dv3;
        points.push_back(new Point(dv1, dv2, dv3));
    }

    int iv1, iv2, iv3;
    for (int i = 0; i < n_face; ++i) {
        infile >> int_temp >> iv1 >> iv2 >> iv3;
        add_face(new Face(iv1, iv2, iv3 ));
    }
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
	points.resize(num_vertices);
	for(int i=0;i<num_vertices;i++){
		Point *p = new Point();
		for(int j=0;j<3;j++){
			if(!parse_float(data, offset, len, p->v[j])){
				return false;
			}
		}
		points[i] = p;
	}
	if(num_faces==0){
		return true;
	}
	faces.resize(num_faces);
	for(int i=0;i<num_faces;i++){
		faces[i] = new Face();
		int fnum = 0;
		if(!parse_int(data, offset, len, fnum) || fnum<3){
			log("invalid face numbers %d", fnum);
			return false;
		}
		faces[i]->vertices.resize(fnum);
		for(int j=0;j<fnum;j++){
			int vindex = 0;
			if(!parse_int(data, offset, len, vindex) || vindex>=points.size()){
				log("invalid vertex index %d", vindex);
				return false;
			}
			faces[i]->vertices[j] = vindex;
		}
	}
	return true;
}

void Polyhedron::add_face(Face *f){
	for(int i=0;i<f->vertices.size();i++){
		pair<int,int> e({f->vertices[i],f->vertices[(i+1)%f->vertices.size()]});
		if(edges.find(e)==edges.end()){
			vector<Face *> fs;
			fs.push_back(f);
			edges[e] = fs;
		}else{
			edges[e].push_back(f);
		}
	}
	faces.push_back(f);
}

void Polyhedron::remove_redundant(){
	map<tuple<int, int, int>, Face *> face_map;
	for(Face *f:faces){

		tuple<int, int, int> tp(f->vertices[0],f->vertices[1],f->vertices[2]);

		if(face_map.find(tp)!=face_map.end()){
			cout<<"redundant face identified"<<endl;
			face_map[tp]->print();
			f->print();
		}else{
			face_map[tp] = f;
		}
	}
}

void Polyhedron::merge_vertex(){
	map<tuple<double, double, double>, int> vertex_map;
	vector<int> vid_map;
	vid_map.resize(points.size());
	int index = 0;
	int removed = 0;
	for(int i=0;i<points.size();i++){
		Point *p = points[i];
		tuple<double, double, double> tp({p->v[0], p->v[1], p->v[2]});
		if(vertex_map.find(tp)==vertex_map.end()){
			vertex_map[tp] = i;
			vid_map[i] = i;
		}else{
			vid_map[i] = vertex_map[tp];
		}
	}
	for(Face *f:faces){
		for(int i=0;i<f->vertices.size();i++){
			f->vertices[i] = vid_map[f->vertices[i]];
		}
	}
}


}


