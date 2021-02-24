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
using namespace std;




class Point
{
public:
    Point(double x1, double x2, double x3) : p1(x1), p2(x2), p3(x3){}
    double p1;   
    double p2;  
    double p3;
};


class Face {
public:
    Face(int p1, int p2, int p3):_idx1(p1), _idx2(p2), _idx3(p3) {

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
    int _idx1, _idx2, _idx3;


};

string check(vector<Face> faces) {
    vector<pair<int, int>> edges;
    for (auto& each : faces) {
        edges.push_back({ each._idx1, each._idx2 });
        edges.push_back({ each._idx2, each._idx3 });
        edges.push_back({ each._idx3, each._idx1 });
    }

    map<pair<int, int>, int> edge_fea;
    for (auto& edge : edges) {
        ++edge_fea[edge];
    }

    vector<pair<int, int>> problem_edges0, problem_edges1;
    for (auto& each : edge_fea) {
        auto rev_k = pair<int, int>{ each.first.second, each.first.first };
        if (!edge_fea.count(rev_k)) {
            problem_edges0.push_back(each.first);
        }
        if (each.second > 1) {
            problem_edges1.push_back(each.first);
        }
    }
    string msg;
    if (!problem_edges0.empty()) {
        msg = "P0 Error";
    }
    else if (!problem_edges1.empty()) {
        msg = "P1 Error ";
    }
    else {
        msg = "";
    }
    set<int> V;
    for (auto& each : faces) {
        V.insert(each._idx1);
        V.insert(each._idx2);
        V.insert(each._idx3);
    }
    int Vs = V.size();
    int F = faces.size();
    int E = edges.size() / 2.0;
    //cout << "V: " << Vs << " F: " << F << " E: " << E << endl;
    //cout << "V+F-E: " << Vs + F - E << endl;
    /*V = len(np.unique(faces))
        F = len(faces)
        E = len(edges) / 2.*/
    //cout << msg << endl;
    return msg;

}

//
// read vertices and faces from original OFF file
//
pair < vector<Point>, vector<Face>> readData(string fp) {
    ifstream infile(fp);
    int n_p, n_face, int_temp;
    string line;

    getline(infile, line);  // OFF
    infile >> n_p >> n_face >> int_temp;
    getline(infile, line); // empty_line
    vector<Point> points;
    double dv1, dv2, dv3;
    for (int i = 0; i < n_p; ++i) {
        infile >> dv1 >> dv2 >> dv3;
        points.emplace_back(Point(dv1, dv2, dv3 ));
    }

    vector<Face> faces;
    int iv1, iv2, iv3;
    for (int i = 0; i < n_face; ++i) {
        infile >> int_temp >> iv1 >> iv2 >> iv3;

        faces.emplace_back(Face(iv1, iv2, iv3 ));
    }
    return { points, faces };
}

pair<vector<Face>, set<int>> depart(vector<Face> faces) {
    vector<Face> inde_poly;
    set<int> ver_poly;


    ver_poly.insert(faces.front()._idx1);
    ver_poly.insert(faces.front()._idx2);
    ver_poly.insert(faces.front()._idx3);


    for (auto each : faces) {
        if (ver_poly.count(each._idx1) || ver_poly.count(each._idx2) || ver_poly.count(each._idx3)) {
            inde_poly.push_back(each);
            ver_poly.insert(each._idx1);
            ver_poly.insert(each._idx2);
            ver_poly.insert(each._idx3);
        }
    }
    return make_pair( inde_poly, ver_poly );
}

bool face_in_pointset(const Face f, const set<int> points) {
    // return True if there is at least one point 
    return points.count(f._idx1) || points.count(f._idx2) || points.count(f._idx3);
}

bool findFaceInSet(Face& ths, vector<Face>& faceSet) {
    return false;
}

vector<Face> simple_copy_if(vector<Face> ths, vector<Face> faceVec) {
    vector<Face> ret;
    for (auto each : ths) {
        /*if (binary_search(faceVec.begin(), faceVec.end(), each)) {
            continue;
        }*/
        if (find(faceVec.begin(), faceVec.end(), each) != faceVec.end()) {
            continue;
        }
        ret.emplace_back(each);
        //if (!face_in_pointset(each, pointSet)) {
        //    ret.emplace_back(each);
        //}
    }
    return ret;
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
// combine polyhedra that connect
vector<int> combine(vector<set<int>> multi_ver_poly, vector<int> unused_idx) {
    
    auto firstUnsedIdx = unused_idx.front();
    unused_idx.erase(unused_idx.begin());

    auto compound = set<int>(multi_ver_poly[firstUnsedIdx].begin(),multi_ver_poly[firstUnsedIdx].end());

    //set<int> compound_idx{ firstUnsedIdx };

    vector<int> compound_idx;
    compound_idx.push_back(firstUnsedIdx);
    int  tmp_num = 0;
   
        for (auto& each : unused_idx) {
            if (each == firstUnsedIdx) {
                continue;
            }
            tmp_num = compound_idx.size();
            //rhttp://www.cplusplus.com/reference/algorithm/set_intersection/
            vector<int> insersectionRet;
            set_intersection(compound.begin(), compound.end(),
                multi_ver_poly[each].begin(), multi_ver_poly[each].end(),
                back_inserter(insersectionRet));

            if (!insersectionRet.empty()) {
                compound.insert(multi_ver_poly[each].begin(), multi_ver_poly[each].end());
                compound_idx.push_back(each);
            }
        }
        //unused_idx = filter_not_in(unused_idx, compound_idx);
    
    return compound_idx;
}


//pair<vector<pair<int, int>>, vector<int>> recognize_ring(vector<pair<int, int>>& problem_edges,
//    vector<int> & unused_idx) {
//    vector<pair<int, int>> ring;
//    vector<int> ring_idx;
//    auto first_idx = *unused_idx.begin();
//    ring.push_back(problem_edges[first_idx]);
//    ring_idx.push_back(first_idx);
//
//
//    while (ring.back().second != ring.front().first) {
//        int count = 0;
//        for (auto& i : unused_idx) {
//            if (i == first_idx) {
//                continue;
//            }
//            if (problem_edges[i].first == ring.back().second) {
//                ring.push_back(problem_edges[i]);
//                ring_idx.push_back(i);
//            }
//            ++count;
//        }
//        if (count > 2 * unused_idx.size()) {
//            break;
//        }
//    }
//    return make_pair(ring, ring_idx);
//
//}




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


vector<int> filter_not_in_set(vector<int> v, set<int> s) {
    vector<int> ret;
    for (auto& each : v) {
        if (s.count(each)) {
            continue;
        }
        ret.push_back(each);
    }
    return ret;
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



int main(int argc, char **argv)
{   
    //####################################################
    //## read vertices and faces from original OFF file ##
    //####################################################
    auto data = readData(argv[1]);
    auto points = data.first;
    auto faces = data.second;

    cout << "Num of Points: " << points.size() << "  Num of Faces: " << faces.size() << endl;


    cout << "Separate a single polyhedron into multiple polyhedra...It may take a while..." << endl;

    //##########################################################
    //## Separate a single polyhedron into multiple polyhedra ##
    //##########################################################


    vector<vector<Face>> multi_poly;  // a list of polyhedra, each item is a list of faces that constitutes a polyhedron
    vector<set<int>> multi_ver_poly;  // a list of points that belong to one polyhedron, each item is an index of a point of a polyhedron
    vector<int> num;
    while (!faces.empty()) {
        auto [inde_poly, ver_poly] = depart(faces); // support by c++17

        //sort(inde_poly.begin(), inde_poly.end());
        auto ret = simple_copy_if(faces, inde_poly);
        faces.clear();
        faces.swap(ret);

        multi_poly.emplace_back(inde_poly);
        multi_ver_poly.push_back(ver_poly);

        num.emplace_back(inde_poly.size());

        //cout << faces.size() << endl;
    }

    cout << "num of used faces:" << accumulate(num.begin(), num.end(),
        decltype(num)::value_type(0)) << endl;;
    cout << "num of initial polyhedra: " << multi_poly.size() << endl;


    // Combine polyhedra that connect

    set<int> used_idx;  // 
    vector<vector<int>> total_poly_idx;

    //set<int> unused_idx; //
    vector<int> unused_idx;
    // init total idx
    vector<int> total_idx(multi_ver_poly.size());
    iota(total_idx.begin(), total_idx.end(), 0);
    while (used_idx.size() != multi_ver_poly.size()) {
        unused_idx.clear();
        for (auto& each : total_idx) {
            if (used_idx.count(each)) {
                continue;
            }
            unused_idx.push_back(each);
        }

        auto compound_idx = combine(multi_ver_poly, unused_idx);
        used_idx.insert(compound_idx.begin(), compound_idx.end());
        total_poly_idx.push_back(compound_idx);
    }

    vector<vector<Face>> total_inde_poly;
    vector<Face> inde_poly;
    for (auto& poly_idx : total_poly_idx) {
        inde_poly.clear();
        for (auto& ieach : poly_idx) {
            for (auto& face : multi_poly[ieach]) {
                inde_poly.push_back(face);
            }
        }
        total_inde_poly.push_back(inde_poly);
    }
    

    //Till now, each polyhedra is independent, and has no connection
    cout << "total num of independent polyhedra: " << total_inde_poly.size() << endl;


    //##############################################################################
    //## process each polyhedron respectively
    //## first recognize rings, and fill holes by adding new faces
    //## then remove redundant faces for those edges which belong to more than 2 faces
    //## finally check whether edges and vertices obey corresponding rules
    //##############################################################################



    vector<Face> fixed_polyhedra;
    for (int i = 0; i < 190; ++i) {
        auto inde_poly = total_inde_poly[i];
        if (inde_poly.size() >= 4) {
            cout << "############## Start processing the" << i << "-th polyhedron ###################################" << endl;
            cout << "original num of inde_poly: " << inde_poly.size() << endl;

            auto inde_poly_ret = fill_holes(inde_poly);

            cout << "total num of inde_poly after filling holes: " << inde_poly_ret.size() << endl;

            auto inde_poly_rev = remove_redundant(inde_poly_ret);

            cout << "total face after remove_redundant: " << inde_poly_rev.size() << endl;
            check(inde_poly_rev);
            copy(inde_poly_rev.begin(), inde_poly_rev.end(), back_inserter(fixed_polyhedra));

            auto msg = check(inde_poly_rev);
            if (msg.size() > 0) {
                cout << "ERROR" << endl;
            }
            else {
                cout << "All fixed." << endl;
            }
        }
    }
    
    cout << "final total num of faces: " << fixed_polyhedra.size() << endl;


    //#####################################
    //## write fixed faces into OFF file ##
    //#####################################

    ofstream of("modi_648518346349489985.OFF");

    of << "OFF" << endl;
    of << points.size() << " " << fixed_polyhedra.size() << "\n" << endl;
    for (auto& each : points) {
        of << each.p1 << " " << each.p2 << " " << each.p3 << endl;
    }
    for (auto& each : fixed_polyhedra) {
        of << "3\t" << each._idx1 << " " << each._idx2 << " " << each._idx3 << endl;
    }

    
  
}
