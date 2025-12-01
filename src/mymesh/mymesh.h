/*
 * mymesh.h
 *
 *  Created on: Nov 2, 2022
 *      Author: teng
 */

#ifndef SRC_MYMESH_MYMESH_H_
#define SRC_MYMESH_MYMESH_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <stack>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "util.h"

using namespace std;
using namespace tdbase;

namespace tmesh {

class Point {
  public:
    float v[3];
    Point(float x1, float x2, float x3) {
        v[0] = x1;
        v[1] = x2;
        v[2] = x3;
    }
    Point(Point* pt) {
        assert(pt);
        for (int i = 0; i < 3; i++) {
            v[i] = pt->v[i];
        }
    };
    Point(){};
};

class Vertex;
class Face;
class Half_Edge;

class Vertex : public Point {
  public:
    Vertex() {}
    Vertex(float v1, float v2, float v3) : Point(v1, v2, v3) {}
    unordered_set<Half_Edge*> half_edges;
    unordered_set<Half_Edge*> opposite_half_edges;
    int id = 0;
    bool added = false;
    bool removable = true;
    int degree() {
        return half_edges.size();
    }
    void print() {
        printf("%f %f %f\n", v[0], v[1], v[2]);
    }
    bool is_removable() {
        return degree() > 2 && removable;
    }
    void setId(int i) {
        id = i;
    }
    int getId() {
        return id;
    }
};

class Half_Edge {
  public:
    Vertex* vertex = NULL;
    Vertex* end_vertex = NULL;
    Face* face = NULL;
    Half_Edge* next = NULL;
    Half_Edge* opposite = NULL;
    Half_Edge(Vertex* v1, Vertex* v2);
    ~Half_Edge();
};

class Face {
  public:
    int id = 0;
    bool added = false;
    vector<Vertex*> vertices;
    unordered_set<Half_Edge*> half_edges;

    Vertex* removedVertexPos = NULL;
    vector<Vertex*> impact_points;
    float conservative_distance = 0.0;
    float progressive_distance = 0.0;
    inline pair<float, float> getHausdorfDistance() {
        return pair<float, float>(conservative_distance, progressive_distance);
    }

  public:
    Face(){};
    ~Face() {
        for (Half_Edge* h : half_edges) {
            delete h;
        }
        half_edges.clear();
        vertices.clear();
    }

    Face(Vertex* v1, Vertex* v2, Vertex* v3) {
        vertices.push_back(v1);
        vertices.push_back(v2);
        vertices.push_back(v3);
    }

    Face(vector<Vertex*>& vs) {
        Half_Edge* prev = NULL;
        Half_Edge* head = NULL;
        for (int i = 0; i < vs.size(); i++) {
            vertices.push_back(vs[i]);
            Vertex* nextv = vs[(i + 1) % vs.size()];
            Half_Edge* hf = new Half_Edge(vs[i], nextv);
            half_edges.insert(hf);
            hf->face = this;
            if (prev != NULL) {
                prev->next = hf;
            } else {
                head = hf;
            }
            if (i == vs.size() - 1) {
                hf->next = head;
            }
            prev = hf;
        }
    }

    void print() {
        printf("totally %ld vertices:\n", vertices.size());
        int idx = 0;
        for (Vertex* v : vertices) {
            printf("%d:\t", idx++);
            v->print();
        }
    }

    void print_off() {
        printf("OFF\n%ld 1 0\n", vertices.size());
        for (Vertex* v : vertices) {
            v->print();
        }
        printf("%ld\t", vertices.size());
        for (int i = 0; i < vertices.size(); i++) {
            printf("%d ", i);
        }
        printf("\n");
    }

    bool equal(const Face& rhs) const {
        if (vertices.size() != rhs.vertices.size()) {
            return false;
        }
        for (int i = 0; i < vertices.size(); i++) {
            if (vertices[i] != rhs.vertices[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator==(const Face& rhs) const {
        return this->equal(rhs);
    }

    int facet_degree() {
        return vertices.size();
    }
    // split the face and make sure the one without v as the new
    Face* split(Vertex* v);
    void remove(Half_Edge* h);
};

class TMesh {
  public:
    int id = 0;
    unordered_set<Vertex*> vertices;
    unordered_set<Face*> faces;
    float low[3];
    float high[3];

    bool owned_data = false;
    char* p_data = NULL;
    size_t dataOffset = 0;

    Vertex* vh_departureConquest[2];

    int i_nbDecimations = 0;
    int i_curDecimationId = 0;

  public:
    TMesh(int i = 0) {
        id = i;
    }
    ~TMesh();

    // I/O
    void load(char* data, bool owned = false);
    void load(string path);
    bool parse(string str);
    bool parse(const char*, size_t);
    void dumpto(string path);
    void print();
    string to_string();
    Vertex* get_vertex(int vseq = 0);

    // element operating
    Face* add_face(vector<Vertex*>& vs);
    Face* remove_vertex(Vertex* v);
    Half_Edge* merge_edges(Vertex* v);

    // compression/decompression
    void reset_states();
    void compress();

    // mesh fixing
    int remove_orphan_vertices();
    void merge_vertex();
    bool fill_holes();
    void remove_redundant();
    vector<TMesh*> depart();
    void evaluate();

    /*
     * statistics
     *
     * */

    size_t size_of_vertices() {
        return vertices.size();
    }

    size_t size_of_facets() {
        return faces.size();
    }

    /**
     * operations to the data buffer, read/write
     */
    void writeCompressedData();
    void readCompressedData();

    void writeBits(uint32_t data, unsigned i_nbBits, char* p_dest, unsigned& i_bitOffset, size_t& offset);
    uint32_t readBits(unsigned i_nbBits, char* p_src, unsigned& i_bitOffset, size_t& offset);
    void writeFloat(float f);
    float readFloat();
    int readInt();
    void writeInt(int i);
    int16_t readInt16();
    void writeInt16(int16_t i);
    uint16_t readuInt16();
    void writeuInt16(uint16_t i);
    unsigned char readChar();
    void writeChar(unsigned char i);

    void writeBaseMesh();
    void readBaseMesh();
    void writeMeshOff(const char psz_filePath[]);
    void writeCurrentOperationMesh(std::string pathPrefix, unsigned i_id);

    vector<pair<float, float>> globalHausdorfDistance;
    pair<float, float> getHausdorfDistance();
    pair<float, float> getNextHausdorfDistance();
};

}  // namespace tmesh
#endif /* SRC_MYMESH_MYMESH_H_ */
