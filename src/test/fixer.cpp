#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <set>
#include <map>
#include <stack>
#include <stdlib.h>
#include "../util/util.h"
using namespace std;


class Point{
public:
    double p1;   
    double p2;  
    double p3;
    Point(double x1, double x2, double x3) : p1(x1), p2(x2), p3(x3){}
    Point(Point *p):p1(p->p1),p2(p->p2),p3(p->p3){};

};

class Face;

class Half_Edge{
public:
	int p1;
	int p2;
	Face *face = NULL;
};

class Face{
public:
	Half_Edge edges[3];
	vector<int> vertices;

	Face(vector<int> &vs){
		vertices.insert(vertices.begin(), vs.begin(), vs.end());
	}
	Face(){};

    Face(int v1, int v2, int v3){

    	vertices.push_back(v1);
    	vertices.push_back(v2);
    	vertices.push_back(v3);

    	edges[0].face = this;
    	edges[0].p1 = v1;
    	edges[0].p2 = v2;

    	edges[1].face = this;
    	edges[1].p1 = v2;
    	edges[1].p2 = v3;

    	edges[2].face = this;
    	edges[2].p1 = v3;
    	edges[2].p2 = v1;
    }

    void print(){
    	for(int v:vertices){
    		printf("%d ",v);
    	}
    	printf("\n");
    }


    bool equal(const Face& rhs) const {
    	if(vertices.size()!=rhs.vertices.size()){
    		return false;
    	}
    	for(int i=0;i<vertices.size();i++){
    		if(vertices[i]!=rhs.vertices[i]){
    			return false;
    		}
    	}
    	return true;
    }
    bool operator==(const Face& rhs) const
    {
        return this->equal(rhs);
    }


};



class Polyhedron{

public:
	int id = 0;
	vector<Point *> points;
	vector<Face *>	faces;
	map<pair<int,int>, vector<Face *>> edges;
	Polyhedron(int i=0){id=i;}
	~Polyhedron();
	void load(string path);
	void dumpto(string path);
	void get_edge();
	void add_face(Face *f);
	void evaluate();
	void merge_vertex();
	void fill_holes();
	void remove_redundant();
	vector<Polyhedron *> depart();
	void print();

};


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
		cout << p->p1 << " " << p->p2 << " " << p->p3 << endl;
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
		of << p->p1 << " " << p->p2 << " " << p->p3 << endl;
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


// to get the largest connected component
vector<Polyhedron *> Polyhedron::depart(){

	vector<vector<int>> connections;
	connections.resize(points.size());
	for(Face *f:faces){
		for(int i=0;i<f->vertices.size();i++){
			connections[f->vertices[i]].push_back(f->vertices[(i+1)%f->vertices.size()]);
			connections[f->vertices[(i+1)%f->vertices.size()]].push_back(f->vertices[i]);
		}
	}
	vector<bool> added;
	for(int i=0;i<points.size();i++){
		added.push_back(false);
	}
	vector<map<int,bool>> ccg;
	while(true){
		map<int,bool> cur;
		int seed = 0;
		for(;seed<points.size();seed++){
			if(!added[seed]){
				break;
			}
		}
		// all added
		if(seed==points.size()){
			break;
		}
		stack<int> ss;
		ss.push(seed);
		while(!ss.empty()){
			seed = ss.top();
			ss.pop();
			cur[seed] = true;
			added[seed] = true;
			// connect all the vertices connect to seed but not visited
			for(int p:connections[seed]){
				if(cur.find(p)==cur.end()){
					ss.push(p);
				}
			}
		}
		ccg.push_back(cur);
	}

	for(vector<int> &con:connections){
		con.clear();
	}
	connections.clear();

	vector<Polyhedron *> polys;
	for(int i=0;i<ccg.size();i++){
		polys.push_back(new Polyhedron(i));
	}
	polys.resize(ccg.size());
	vector<int> pid_map;
	vector<int> cg_map;
	pid_map.resize(points.size());
	cg_map.resize(points.size());

	for(int i=0;i<ccg.size();i++){
		for(auto a:ccg[i]){
			assert(a.first<points.size());
			pid_map[a.first] = polys[i]->points.size();
			polys[i]->points.push_back(new Point(points[a.first]));
			cg_map[a.first] = i;
		}
		ccg[i].clear();
	}
	ccg.clear();

	// assign all the faces
	for(Face *f:faces){
		Face *nf = new Face();
		int cg_id = cg_map[f->vertices[0]];
		for(int v:f->vertices){
			assert(cg_map[v]==cg_id);
			nf->vertices.push_back(pid_map[v]);
		}
		polys[cg_id]->add_face(nf);
	}

	pid_map.clear();
	cg_map.clear();
	return polys;
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
		tuple<double, double, double> tp({p->p1, p->p2, p->p3});
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

void Polyhedron::fill_holes(){

	// identify the single edges
	vector<pair<int,int>> single_edges;
	for(auto &entry:edges){
		pair<int,int> rev({entry.first.second, entry.first.first});
		if(edges.find(rev)==edges.end()){
			single_edges.push_back(entry.first);
		}
	}
	// no single_edges
	if(single_edges.size()==0){
		return;
	}

	// do another connect component graph find
	map<int,bool> assigned;
	map<int,int> connect;
	for(pair<int, int> &e:single_edges){
		// should be connected in reverse order
		connect[e.second] = e.first;
		assigned[e.first] = false;
		assigned[e.second] = false;
	}
	cout<<"boundary edge: "<<single_edges.size()<<endl;
	cout<<"boundary vertex: "<<assigned.size()<<endl;

	while(true){
		int curv = -1;
		// find a vertex
		for(auto &v:assigned){
			if(!v.second){
				curv = v.first;
				break;
			}
		}
		// all assigned
		if(curv==-1){
			break;
		}
		Face *f = new Face();
		bool valid = true;
		while(true){
			f->vertices.push_back(curv);
			assigned[curv] = true;
			//cout<<curv<<endl;
			// cannot find a loop in this round
			if(connect.find(curv)==connect.end()){
				valid = false;
				break;
			}
			// forward to next
			curv = connect[curv];
			// should form a loop now, current vertex must be the first
			if(assigned[curv]){
				valid = (curv == f->vertices[0]);
				//assert(curv == f->vertices[0]);
				break;
			}
		}
		if(valid){
			add_face(f);
		}else{
			delete f;
		}
	}
}

void Polyhedron::evaluate(){
    int V = points.size();
    int F = faces.size();
    int E = 0;
    for(Face *f:faces){
    	E += f->vertices.size();
    }
    E /= 2;
    cout<<"Polyhedron: "<<id<<endl;
    cout << "V: " << V << " F: " << F << " E: " << E << endl;
    cout << "V+F-E: " << V + F - E << endl;
}



int main(int argc, char **argv){
	Polyhedron *poly = new Polyhedron();
	poly->load(argv[1]);
	vector<Polyhedron *> ind = poly->depart();

	char path[256];
	for(int i=0;i<ind.size();i++){
		if(ind[i]->points.size()>100){
			ind[i]->evaluate();
			//ind[i]->remove_redundant();
			ind[i]->merge_vertex();
			ind[i]->fill_holes();
			ind[i]->evaluate();
			//ind[i]->fill_holes();
			sprintf(path,"/gisdata/updated_%ld_%ld.OFF",ind[i]->points.size(),ind[i]->faces.size());
			ind[i]->dumpto(path);
		}
		delete ind[i];
	}

	ind.clear();
	delete poly;
}

//
//
//
//vector<int> filter_not_in(vector<int> v, vector<int> s) {
//    vector<int> ret;
//    for (auto& each : v) {
//        if (find(s.begin(), s.end(), each) != s.end()) {
//            continue;
//        }
//        ret.push_back(each);
//    }
//    return ret;
//}
//
//pair<vector<pair<int, int>>, vector<int>> recognize_ring(vector<pair<int, int>> problem_edges,
//    vector<int> unused_idx) {
//
//    vector<pair<int, int>> ring;
//    vector<int> ring_idx;
//    auto first_idx = *unused_idx.begin();
//    ring.push_back(problem_edges[first_idx]);
//    ring_idx.push_back(first_idx);
//
//    unused_idx.erase(unused_idx.begin());
//
//
//    while (ring.back().second != ring.front().first) {
//        int count = 0;
//        for (auto& i : unused_idx) {
//           /* if (i == first_idx) {
//                continue;
//            }*/
//            if (problem_edges[i].first == ring.back().second) {
//                ring.push_back(problem_edges[i]);
//                ring_idx.push_back(i);
//            }
//
//        }
//        auto ret = filter_not_in(unused_idx, ring_idx);
//        unused_idx.clear();
//        unused_idx.swap(ret);
//
//    }
//    return make_pair(ring, ring_idx);
//}
//
//
//
//
//vector<Face> add_face(vector<pair<int, int>> ring) {
//    vector<Face> new_face;
//    pair<int, int> tmp = make_pair(ring.front().second, ring.front().first);
//    while (new_face.size() != ring.size() - 2) {
//        for (auto& edge : ring) {
//            if (edge.second == tmp.second) {
//                //new_face.push_back({ tmp.first, tmp.second, edge.first });
//                new_face.push_back(Face(tmp.first, tmp.second, edge.first));
//                tmp = make_pair(ring.front().second, edge.first);
//            }
//        }
//    }
//    return new_face;
//}
//
//vector<Face> fill_holes(vector<Face> faces) {
//    vector<pair<int, int>> edges;
//    for (auto& e : faces) {
//        edges.push_back(make_pair( e.idx[0], e.idx[1] ));
//        edges.push_back(make_pair( e.idx[1], e.idx[2] ));
//        edges.push_back(make_pair( e.idx[2], e.idx[0] ));
//    }
//    map<pair<int, int>, int> edge_fea;
//    for (auto& e : edges) {
//        ++edge_fea[e];
//    }
//
//    vector<pair<int, int>> problem_edges;
//    for (auto& e : edge_fea) {
//        pair<int, int> rever_edge = make_pair(e.first.second, e.first.first);
//        if (edge_fea.count(rever_edge) == 0) {
//            problem_edges.emplace_back(e.first);
//        }
//    }
//
//    vector< vector<pair<int, int>>> total_ring;
//    vector<int> total_idx(problem_edges.size());
//    iota(total_idx.begin(), total_idx.end(), 0);
//    vector<int> used_idx;
//
//    while (used_idx.size() != problem_edges.size()) {
//        auto unused_idx = filter_not_in(total_idx, used_idx);
//
//        //set<int> unused_idx_set(unused_idx.begin(), unused_idx.end());
//        auto [ring, ring_idx] = recognize_ring(problem_edges, unused_idx);
//
//        copy(ring_idx.begin(), ring_idx.end(), back_inserter(used_idx));
//        total_ring.push_back(ring);
//    }
//
//    vector<Face> add_faces;
//    for (auto& ring : total_ring) {
//        auto new_face = add_face(ring);
//        copy(new_face.begin(), new_face.end(), back_inserter(add_faces));
//    }
//    //cout << "num of added faces: " << add_faces.size() << endl;
//
//
//    vector<Face> face_ret(faces.begin(), faces.end());
//    copy(add_faces.begin(), add_faces.end(), back_inserter(face_ret));
//
//    //cout << "check after add face " << endl;
//    //check(face_ret);
//    return face_ret;
//}
//
//
//

//vector<Face> remove_redundant( vector<Face> faces) {
//    vector<pair<int, int>> edges;
//    //sort(faces.begin(), faces.end());
//    for (auto& face : faces) {
//        edges.push_back(make_pair( face.idx[0], face.idx[1] ));
//        edges.push_back(make_pair( face.idx[1], face.idx[2] ));
//        edges.push_back(make_pair( face.idx[2], face.idx[0] ));
//    }
//
//    map<pair<int, int>, int> edge_fea;
//    for (auto& edge : edges) {
//        ++edge_fea[edge];
//    }
//
//    vector<pair<int, int>> problem_edges;
//    for (auto& e : edge_fea) {
//        if (e.second > 1) {
//            for (int i = 0; i < e.second - 1; ++i) {
//                problem_edges.push_back(e.first);
//            }
//        }
//    }
//
//    //cout << "num of edges that belong to more than two faces: " << problem_edges.size() << endl;
//
//    vector<vector<pair<int, int>>> pair_prob_edges;
//    for (auto& edge : problem_edges) {
//        pair<int, int> rev_edge{ edge.second, edge.first };
//        if (find(problem_edges.begin(), problem_edges.end(), rev_edge) != problem_edges.end()) {
//
//            auto insert_item = vector<pair<int, int>>{ edge, rev_edge };
//            if (insert_item[0].first > insert_item[1].first) {
//                swap(insert_item[0], insert_item[1]);
//            }
//            else if (insert_item[0].second > insert_item[1].second) {
//                swap(insert_item[0], insert_item[1]);
//            }
//            //sort(insert_item.begin(), insert_item.end(), [](const pair<int, int>& a, pair<int, int>& b) {a.first < b.first ? true : a.second < b.second; });
//            pair_prob_edges.push_back(insert_item);
//        }
//    }
//
//    vector<vector<pair<int, int>>> pair_prob_edges_;
//    for (auto& each : pair_prob_edges) {
//         // find if not in
//        if (find(pair_prob_edges_.begin(), pair_prob_edges_.end(), each) == pair_prob_edges_.end()) {
//            pair_prob_edges_.push_back(each);
//        }
//    }
//
//    //cout << "num of paired edges: " << pair_prob_edges_.size() << endl;
//
//    vector<vector<pair<int, int>>> face_edges;
//    for (auto& face : faces) {
//        face_edges.push_back(
//            { make_pair(face.idx[0], face.idx[1]),
//            make_pair(face.idx[1], face.idx[2]),
//            make_pair(face.idx[2], face.idx[0]), }
//        );
//    }
//
//    vector<Face> redundant;
//    for (auto& pair_edge : pair_prob_edges_) {
//        vector<Face> focused;
//        for (int i = 0; i < face_edges.size(); ++i) {
//            auto face_edge = face_edges[i];
//
//            if (
//                find(face_edge.begin(), face_edge.end(), pair_edge[0]) != face_edge.end() ||
//                find(face_edge.begin(), face_edge.end(), pair_edge[1]) != face_edge.end()
//                ) {
//                focused.push_back(faces[i]);
//            }
//        }
//        for (auto& face : focused) {
//            auto face1 = Face{ face.idx[2], face.idx[1], face.idx[0] };
//            if (
//                find(focused.begin(), focused.end(),face1) != focused.end()
//                ) {
//                redundant.push_back(face);
//                redundant.push_back(face1);
//                break;
//            }
//
//            auto face2 = Face{ face.idx[0], face.idx[2], face.idx[1] };
//            if (
//                find(focused.begin(), focused.end(), face2) != focused.end()
//                ) {
//                redundant.push_back(face);
//                redundant.push_back(face2);
//                break;
//            }
//
//            auto face3 = Face{ face.idx[1], face.idx[0], face.idx[2] };
//            if (
//                find(focused.begin(), focused.end(), face3) != focused.end()
//                ) {
//                redundant.push_back(face);
//                redundant.push_back(face3);
//                break;
//            }
//
//        }
//
//    }
//
//	vector<Face> faces_ret;
//	for (auto& e : faces) {
//		// ignore found item
//		if (find(redundant.begin(), redundant.end(), e) != redundant.end()) {
//			continue;
//		}
//		faces_ret.push_back(e);
//	}
//
//    return faces_ret;
//
//}
//
//
//
//int main1(int argc, char **argv){
//    //####################################################
//    //## read vertices and faces from original OFF file ##
//    //####################################################
//    vector<Point> points;
//    vector<Face> faces;
//
//    vector<vector<Face>> total_inde_poly;
//
//    //##############################################################################
//    //## process each polyhedron respectively
//    //## first recognize rings, and fill holes by adding new faces
//    //## then remove redundant faces for those edges which belong to more than 2 faces
//    //## finally check whether edges and vertices obey corresponding rules
//    //#############################################################################
//
//
//    vector<Face> fixed_polyhedra;
//    for (int i = 0; i < total_inde_poly.size(); ++i) {
//        auto inde_poly = total_inde_poly[i];
//        if (inde_poly.size() >= 4) {
//            cout << "############## Start processing the" << i << "-th polyhedron ###################################" << endl;
//            cout << "original num of inde_poly: " << inde_poly.size() << endl;
//
//            auto inde_poly_ret = fill_holes(inde_poly);
//
//            cout << "total num of inde_poly after filling holes: " << inde_poly_ret.size() << endl;
//
//            auto inde_poly_rev = remove_redundant(inde_poly_ret);
//
//        }
//    }
//    return 0;
//}
