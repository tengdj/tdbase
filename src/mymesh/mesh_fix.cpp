/*
 * mesh_fix.cpp
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */


#include "mymesh.h"


namespace tmesh{

// to get the largest connected component
vector<TMesh *> TMesh::depart(){
	vector<TMesh *> polys;

	vector<bool> added;
	for(int i=0;i<vertices.size();i++){
		added.push_back(false);
	}

	vector<map<Vertex *, bool>> ccg;
	while(true){
		map<Vertex *,bool> cur;
		Vertex *seed;
		for(int i=0;i<vertices.size();i++){
			if(!added[i]){
				seed = vertices[i];
				break;
			}
		}

		// all added
		if(seed == NULL){
			break;
		}
		stack<Vertex *> ss;
		ss.push(seed);
		while(!ss.empty()){
			seed = ss.top();
			ss.pop();
			cur[seed] = true;
			added[seed->id] = true;
			// connect all the vertices connect to seed but not visited
			for(Half_Edge *h:seed->half_edges){
				if(cur.find(h->end_vertex)==cur.end()){
					ss.push(h->end_vertex);
				}
			}
		}
		ccg.push_back(cur);
	}

	for(int i=0;i<ccg.size();i++){
		polys.push_back(new TMesh(i));
	}

	for(int i=0;i<ccg.size();i++){
		for(auto &a:ccg[i]){
			assert(a.first->id<vertices.size());
			a.first->id = polys[i]->vertices.size();
			polys[i]->vertices.emplace(a.first);
		}
		ccg[i].clear();
	}
	ccg.clear();

	// retrieve the faces and half_edges from the vertices
	for(TMesh *p:polys){
		unordered_set<Half_Edge *> hedges;
		unordered_set<Face *> faces;
		for(Vertex *v:p->vertices){
			for(Half_Edge *h:v->half_edges){
				if(faces.find(h->face)==faces.end()){
					faces.insert(h->face);
					p->faces.emplace(h->face);
				}
			}
		}
		hedges.clear();
		faces.clear();
	}

	return polys;
}

bool TMesh::fill_holes(){

//	while(true){
//		// identify the singsingle_edgele edges
//		Half_Edge *single_edge = NULL;
//		for(Half_Edge *h:half_edges){
//			if(h->opposite == NULL){
//				single_edge = h;
//			}
//		}
//
//		// no single_edges
//		if(single_edge == NULL){
//			return true;
//		}
//
//
//
//	}
	return true;

}

void TMesh::remove_redundant(){
//	map<tuple<int, int, int>, Face *> face_map;
//	for(Face *f:faces){
//
//		tuple<int, int, int> tp(f->vertices[0],f->vertices[1],f->vertices[2]);
//
//		if(face_map.find(tp)!=face_map.end()){
//			cout<<"redundant face identified"<<endl;
//			face_map[tp]->print();
//			f->print();
//		}else{
//			face_map[tp] = f;
//		}
//	}
}

void TMesh::merge_vertex(){
//	map<tuple<double, double, double>, int> vertex_map;
//	vector<int> vid_map;
//	vid_map.resize(points.size());
//	int index = 0;
//	int removed = 0;
//	for(int i=0;i<points.size();i++){
//		Point *p = points[i];
//		tuple<double, double, double> tp({p->v[0], p->v[1], p->v[2]});
//		if(vertex_map.find(tp)==vertex_map.end()){
//			vertex_map[tp] = i;
//			vid_map[i] = i;
//		}else{
//			vid_map[i] = vertex_map[tp];
//		}
//	}
//	for(Face *f:faces){
//		for(int i=0;i<f->vertices.size();i++){
//			f->vertices[i] = vid_map[f->vertices[i]];
//		}
//	}
}

void TMesh::evaluate(){

	for(Face *f:faces){
		for(Half_Edge *he:f->half_edges){
			assert(he->opposite && "hole exist");
		}
	}

    int V = vertices.size();
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
