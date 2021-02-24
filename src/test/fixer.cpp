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
    int _idx1, _idx2, _idx3;

    Face(int p1, int p2, int p3):_idx1(p1), _idx2(p2), _idx3(p3) {
    	edges[0].face = this;
    	edges[0].p1 = _idx1;
    	edges[0].p2 = _idx2;

    	edges[1].face = this;
    	edges[1].p1 = _idx2;
    	edges[1].p2 = _idx3;

    	edges[2].face = this;
    	edges[2].p1 = _idx3;
    	edges[2].p2 = _idx1;
    }

    bool operator<(const Face& rhs) const
    {
        if (_idx1 <= rhs._idx1) {
            return true;
        }
        else if (_idx1 > rhs._idx1) {
            return false;
        }
        else if (_idx2 <= rhs._idx2)
        {
            return true;
        }
        else if (_idx2 > rhs._idx2) {
            return false;
        }
        else if (_idx3 <= rhs._idx3) {
            return true;
        }
        else if (_idx3 > rhs._idx3) {
            return false;
        }
        return false;
    }
    bool equal(const Face& rhs) {
        return rhs._idx1 == _idx1&& rhs._idx2 == _idx2&& rhs._idx3 == _idx3;
    }
    bool operator==(const Face& rhs) const
    {
        return rhs._idx1 == _idx1 && rhs._idx2 == _idx2 && rhs._idx3 == _idx3;
    }


};



class Polyhedron{

public:
	int id = 0;
	vector<Point *> points;
	vector<Face *>	faces;
	Polyhedron(int i=0){id=i;}
	void load(string path);
	void dumpto(string path);
	void evaluate();
	vector<Polyhedron *> depart();
	~Polyhedron(){
		for(Point *p:points){
			delete p;
		}
		for(Face *f:faces){
			delete f;
		}
		points.clear();
		faces.clear();
	}
};


void Polyhedron::dumpto(string fp){

	ofstream of(fp.c_str());
	of << "OFF" << endl;
	of << points.size() << " " << faces.size() << " 0\n" << endl;
	for(Point *p:points){
		of << p->p1 << " " << p->p2 << " " << p->p3 << endl;
	}
	for(Face *f:faces){
		of << "3\t" << f->_idx1 << " " << f->_idx2 << " " << f->_idx3 << endl;
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
        faces.push_back(new Face(iv1, iv2, iv3 ));
    }
}


// to get the largest connected component
vector<Polyhedron *> Polyhedron::depart(){

	vector<vector<int>> connections;
	connections.resize(points.size());
	for(Face *f:faces){
		connections[f->_idx1].push_back(f->_idx2);
		connections[f->_idx1].push_back(f->_idx3);
		connections[f->_idx2].push_back(f->_idx1);
		connections[f->_idx2].push_back(f->_idx3);
		connections[f->_idx3].push_back(f->_idx1);
		connections[f->_idx3].push_back(f->_idx2);
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
		int cg_id = cg_map[f->_idx1];
		assert(cg_map[f->_idx1]==cg_map[f->_idx2]&&cg_map[f->_idx1]==cg_map[f->_idx3]);
		polys[cg_id]->faces.push_back(new Face(pid_map[f->_idx1],pid_map[f->_idx2],pid_map[f->_idx3]));
	}

	pid_map.clear();
	cg_map.clear();
	return polys;
}

void Polyhedron::evaluate(){
	map<pair<int, int>, vector<Face *>> edges;

	int num_edges = 0;
	for(Face *f:faces){
		pair<int,int> e({f->_idx1,f->_idx2});
		if(edges.find(e)==edges.end()){
			vector<Face *> fs;
			fs.push_back(f);
			edges[e] = fs;
		}else{
			edges[e].push_back(f);
		}
		e = pair<int,int>({f->_idx2,f->_idx3});
		if(edges.find(e)==edges.end()){
			vector<Face *> fs;
			fs.push_back(f);
			edges[e] = fs;
		}else{
			edges[e].push_back(f);
		}
		e = pair<int,int>({f->_idx3,f->_idx1});
		if(edges.find(e)==edges.end()){
			vector<Face *> fs;
			fs.push_back(f);
			edges[e] = fs;
		}else{
			edges[e].push_back(f);
		}
		num_edges += 3;
	}

	vector<pair<int,int>> pbl1_face;
	vector<pair<int,int>> pbl2_face;

	for(auto &entry:edges){
		bool wrong = false;
		if(entry.second.size()!=1){
			pbl1_face.push_back(entry.first);
		}
		if(edges.find(pair<int, int>({entry.first.second, entry.first.first}))==edges.end()){
			pbl2_face.push_back(entry.first);
		}
	}

	cout<<pbl1_face.size()<<" "<<pbl2_face.size()<<endl;

    int V = points.size();
    int F = faces.size();
    int E = num_edges/2.0;
    cout<<id<<endl;
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
			sprintf(path,"departed_%d.OFF",i);
			ind[i]->evaluate();
			//ind[i]->dumpto(path);
		}
		delete ind[i];
	}

	ind.clear();
	delete poly;
}




vector<int> filter_not_in(vector<int> v, vector<int> s) {
    vector<int> ret;
    for (auto& each : v) {
        if (find(s.begin(), s.end(), each) != s.end()) {
            continue;
        }
        ret.push_back(each);
    }
    return ret;
}

pair<vector<pair<int, int>>, vector<int>> recognize_ring(vector<pair<int, int>> problem_edges,
    vector<int> unused_idx) {

    vector<pair<int, int>> ring;
    vector<int> ring_idx;
    auto first_idx = *unused_idx.begin();
    ring.push_back(problem_edges[first_idx]);
    ring_idx.push_back(first_idx);

    unused_idx.erase(unused_idx.begin());


    while (ring.back().second != ring.front().first) {
        int count = 0;
        for (auto& i : unused_idx) {
           /* if (i == first_idx) {
                continue;
            }*/
            if (problem_edges[i].first == ring.back().second) {
                ring.push_back(problem_edges[i]);
                ring_idx.push_back(i);
            }

        }
        auto ret = filter_not_in(unused_idx, ring_idx);
        unused_idx.clear();
        unused_idx.swap(ret);

    }
    return make_pair(ring, ring_idx);
}

vector<Face> filter_not_in_face(vector<Face> v, vector<Face> s) {
    vector<Face> ret;
    for (auto& e : v) {
        // ignore found item
        if (find(s.begin(), s.end(), e) != s.end()) {
            continue;
        }
        ret.push_back(e);
    }
    return ret;
}


vector<Face> add_face(vector<pair<int, int>> ring) {
    vector<Face> new_face;
    pair<int, int> tmp = make_pair(ring.front().second, ring.front().first);
    while (new_face.size() != ring.size() - 2) {
        for (auto& edge : ring) {
            if (edge.second == tmp.second) {
                //new_face.push_back({ tmp.first, tmp.second, edge.first });
                new_face.push_back(Face(tmp.first, tmp.second, edge.first));
                tmp = make_pair(ring.front().second, edge.first);
            }
        }
    }
    return new_face;
}

vector<Face> fill_holes(vector<Face> faces) {
    vector<pair<int, int>> edges;
    for (auto& e : faces) {
        edges.push_back(make_pair( e._idx1, e._idx2 ));
        edges.push_back(make_pair( e._idx2, e._idx3 ));
        edges.push_back(make_pair( e._idx3, e._idx1 ));
    }
    map<pair<int, int>, int> edge_fea;
    for (auto& e : edges) {
        ++edge_fea[e];
    }

    vector<pair<int, int>> problem_edges;
    for (auto& e : edge_fea) {
        pair<int, int> rever_edge = make_pair(e.first.second, e.first.first);
        if (edge_fea.count(rever_edge) == 0) {
            problem_edges.emplace_back(e.first);
        }
    }

    vector< vector<pair<int, int>>> total_ring;
    vector<int> total_idx(problem_edges.size());
    iota(total_idx.begin(), total_idx.end(), 0);
    vector<int> used_idx;

    while (used_idx.size() != problem_edges.size()) {
        auto unused_idx = filter_not_in(total_idx, used_idx);

        //set<int> unused_idx_set(unused_idx.begin(), unused_idx.end());
        auto [ring, ring_idx] = recognize_ring(problem_edges, unused_idx);

        copy(ring_idx.begin(), ring_idx.end(), back_inserter(used_idx));
        total_ring.push_back(ring);
    }

    vector<Face> add_faces;
    for (auto& ring : total_ring) {
        auto new_face = add_face(ring);
        copy(new_face.begin(), new_face.end(), back_inserter(add_faces));
    }
    //cout << "num of added faces: " << add_faces.size() << endl;


    vector<Face> face_ret(faces.begin(), faces.end());
    copy(add_faces.begin(), add_faces.end(), back_inserter(face_ret));

    //cout << "check after add face " << endl;
    //check(face_ret);
    return face_ret;
}


vector<Face> remove_redundant( vector<Face> faces) {
    vector<pair<int, int>> edges;
    //sort(faces.begin(), faces.end());
    for (auto& face : faces) {
        edges.push_back(make_pair( face._idx1, face._idx2 ));
        edges.push_back(make_pair( face._idx2, face._idx3 ));
        edges.push_back(make_pair( face._idx3, face._idx1 ));
    }

    map<pair<int, int>, int> edge_fea;
    for (auto& edge : edges) {
        ++edge_fea[edge];
    }

    vector<pair<int, int>> problem_edges;
    for (auto& e : edge_fea) {
        if (e.second > 1) {
            for (int i = 0; i < e.second - 1; ++i) {
                problem_edges.push_back(e.first);
            }
        }
    }

    //cout << "num of edges that belong to more than two faces: " << problem_edges.size() << endl;

    vector<vector<pair<int, int>>> pair_prob_edges;
    for (auto& edge : problem_edges) {
        pair<int, int> rev_edge{ edge.second, edge.first };
        if (find(problem_edges.begin(), problem_edges.end(), rev_edge) != problem_edges.end()) {

            auto insert_item = vector<pair<int, int>>{ edge, rev_edge };
            if (insert_item[0].first > insert_item[1].first) {
                swap(insert_item[0], insert_item[1]);
            }
            else if (insert_item[0].second > insert_item[1].second) {
                swap(insert_item[0], insert_item[1]);
            }
            //sort(insert_item.begin(), insert_item.end(), [](const pair<int, int>& a, pair<int, int>& b) {a.first < b.first ? true : a.second < b.second; });
            pair_prob_edges.push_back(insert_item);
        }
    }

    vector<vector<pair<int, int>>> pair_prob_edges_;
    for (auto& each : pair_prob_edges) {
         // find if not in
        if (find(pair_prob_edges_.begin(), pair_prob_edges_.end(), each) == pair_prob_edges_.end()) {
            pair_prob_edges_.push_back(each);
        }
    }

    //cout << "num of paired edges: " << pair_prob_edges_.size() << endl;

    vector<vector<pair<int, int>>> face_edges;
    for (auto& face : faces) {
        face_edges.push_back(
            { make_pair(face._idx1, face._idx2),
            make_pair(face._idx2, face._idx3),
            make_pair(face._idx3, face._idx1), }
        );
    }

    vector<Face> redundant;
    for (auto& pair_edge : pair_prob_edges_) {
        vector<Face> focused;
        for (int i = 0; i < face_edges.size(); ++i) {
            auto face_edge = face_edges[i];

            if (
                find(face_edge.begin(), face_edge.end(), pair_edge[0]) != face_edge.end() ||
                find(face_edge.begin(), face_edge.end(), pair_edge[1]) != face_edge.end()
                ) {
                focused.push_back(faces[i]);
            }
        }
        for (auto& face : focused) {
            auto face1 = Face{ face._idx3, face._idx2, face._idx1 };
            if (
                find(focused.begin(), focused.end(),face1) != focused.end()
                ) {
                redundant.push_back(face);
                redundant.push_back(face1);
                break;
            }

            auto face2 = Face{ face._idx1, face._idx3, face._idx2 };
            if (
                find(focused.begin(), focused.end(), face2) != focused.end()
                ) {
                redundant.push_back(face);
                redundant.push_back(face2);
                break;
            }

            auto face3 = Face{ face._idx2, face._idx1, face._idx3 };
            if (
                find(focused.begin(), focused.end(), face3) != focused.end()
                ) {
                redundant.push_back(face);
                redundant.push_back(face3);
                break;
            }

        }

    }

    auto faces_ret = filter_not_in_face(faces, redundant);
    return faces_ret;

}



int main1(int argc, char **argv){
    //####################################################
    //## read vertices and faces from original OFF file ##
    //####################################################
    vector<Point> points;
    vector<Face> faces;

    vector<vector<Face>> total_inde_poly;

    //##############################################################################
    //## process each polyhedron respectively
    //## first recognize rings, and fill holes by adding new faces
    //## then remove redundant faces for those edges which belong to more than 2 faces
    //## finally check whether edges and vertices obey corresponding rules
    //#############################################################################


    vector<Face> fixed_polyhedra;
    for (int i = 0; i < total_inde_poly.size(); ++i) {
        auto inde_poly = total_inde_poly[i];
        if (inde_poly.size() >= 4) {
            cout << "############## Start processing the" << i << "-th polyhedron ###################################" << endl;
            cout << "original num of inde_poly: " << inde_poly.size() << endl;

            auto inde_poly_ret = fill_holes(inde_poly);

            cout << "total num of inde_poly after filling holes: " << inde_poly_ret.size() << endl;

            auto inde_poly_rev = remove_redundant(inde_poly_ret);

        }
    }
    return 0;
}
