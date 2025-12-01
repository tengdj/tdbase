/*
 * mymesh.cpp
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */

#include "mymesh.h"
#include <queue>
#include <sstream>

namespace tmesh {

TMesh::~TMesh() {
    for (Face* f : faces) {
        delete f;
    }
    for (Vertex* p : vertices) {
        assert(p->half_edges.size() == 0 && p->opposite_half_edges.size() == 0);
        delete p;
    }
    vertices.clear();
    faces.clear();
    if (owned_data) {
        delete p_data;
    }
}

Vertex* TMesh::get_vertex(int vseq) {
    assert(vseq >= 0 && vseq < vertices.size());
    for (Vertex* v : vertices) {
        if (vseq-- == 0) {
            return v;
        }
    }
    assert(false);
    return NULL;
}

Face* TMesh::add_face(vector<Vertex*>& vs) {
    Face* f = new Face(vs);
    faces.insert(f);
    return f;
}

Face* TMesh::remove_vertex(Vertex* v) {
    vector<Face*> face_removed;
    for (Half_Edge* h : v->half_edges) {
        face_removed.push_back(h->face);
    }

    unordered_set<Half_Edge*> boundary;
    for (Face* f : face_removed) {
        for (Half_Edge* h : f->half_edges) {
            if (h->end_vertex != v && h->vertex != v && h->opposite) {
                if (h->opposite) {
                    boundary.insert(h->opposite);
                }
            }
        }
        faces.erase(f);
        delete f;
    }
    face_removed.clear();

    vertices.erase(v);

    assert(v->half_edges.size() == 0 && v->opposite_half_edges.size() == 0);
    vector<Vertex*> ordered_boundary;

    Half_Edge* cur = (*boundary.begin());
    Vertex* head = cur->end_vertex;
    do {
        boundary.erase(cur);
        ordered_boundary.push_back(cur->end_vertex);
        for (auto h : boundary) {
            if (h->end_vertex == cur->vertex) {
                cur = h;
                break;
            }
        }
    } while (!boundary.empty());
    Face* nface = add_face(ordered_boundary);
    nface->added = true;
    ordered_boundary.clear();

    return nface;
}

void TMesh::reset_states() {
    for (Vertex* v : vertices) {
        v->added = false;
        v->removable = true;
    }

    for (Face* f : faces) {
        f->added = false;
    }
}

void TMesh::compress() {
    Vertex* seed = get_vertex(0);
    queue<Vertex*> wq;
    wq.push(seed);
    while (!wq.empty()) {
        Vertex* v = wq.front();
        wq.pop();
        for (Half_Edge* he : v->half_edges) {
            if (v->removable) {
                he->end_vertex->removable = false;
            }
            if (!he->end_vertex->added) {
                he->end_vertex->added = true;
                wq.push(he->end_vertex);
            }
        }
        if (v->is_removable()) {
            vector<Half_Edge*> hes;
            for (Half_Edge* he : v->half_edges) {
                hes.push_back(he);
            }
            assert(hes.size() > 1);
            for (Half_Edge* he : hes) {
                Face* f = he->face;
                assert(f);
                if (f->facet_degree() > 3) {
                    Face* nf = f->split(v);
                    if (nf) {
                        faces.insert(nf);
                    }
                }
            }
            remove_vertex(v);
        }
    }
}

int TMesh::remove_orphan_vertices() {
    vector<Vertex*> orphan;
    for (auto v : vertices) {
        if (v->half_edges.size() == 0) {
            assert(v->opposite_half_edges.size() == 0);
            orphan.push_back(v);
        }
    }
    int ret = orphan.size();
    for (auto v : orphan) {
        vertices.erase(v);
        delete v;
    }
    orphan.clear();
    return ret;
}

}  // namespace tmesh
