/*
 * mesh_fix.cpp
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */


#include "mymesh.h"


namespace hispeed{

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


}
