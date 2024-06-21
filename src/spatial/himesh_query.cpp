/*
 * himesh_query.cpp
 *
 *  Created on: Sep 24, 2022
 *      Author: teng
 */

#include "himesh.h"
#include "geometry.h"

namespace tdbase{

query_context global_ctx;

void HiMesh::print_triangles(float *triangle, size_t size){
	printf("OFF\n%ld %ld 0\n\n",size*3,size);
	for(size_t i=0;i<size*3;i++){
		printf("%f %f %f\n", *(triangle+3*i), *(triangle+3*i+1), *(triangle+3*i+2));
	}
	for(size_t i=0;i<size;i++){
		printf("3\t%ld %ld %ld\n",i*3,i*3+1,i*3+2);
	}
	printf("\n");
}

bool HiMesh::intersect(HiMesh *target){
	float *tri1, *tri2;
	size_t s1 = fill_triangles(tri1);
	size_t s2 = target->fill_triangles(tri2);
	result_container res = MeshInt(tri1, tri2, s1, s2);
	if(res.intersected && global_ctx.verbose>=1){
		print_triangles(tri1+res.p1*9, 1);
		print_triangles(tri2+res.p2*9, 1);
	}
	delete []tri1;
	delete []tri2;
	return res.intersected;
}

float HiMesh::distance(HiMesh *target){
	result_container res;
	float *tri1, *tri2;
	size_t s1 = fill_triangles(tri1);
	size_t s2 = target->fill_triangles(tri2);
	res = MeshDist(tri1, tri2, s1, s2);
	if(global_ctx.verbose>=1){
		print_triangles(tri1+res.p1*9, 1);
		print_triangles(tri2+res.p2*9, 1);
	}
	delete []tri1;
	delete []tri2;

	return res.distance;
}

bool HiMesh::intersect_tree(HiMesh *target){
	list<Segment> segments = get_segments();
	assert(segments.size()>0);
    for(Segment &s:segments){
    	if(target->get_aabb_tree_triangle()->do_intersect(s)){
    		return true;
    	}
    }
    return false;
}

float HiMesh::distance_tree(const Point &p){
	FT sqd = get_aabb_tree_triangle()->squared_distance(p);
	return sqrt((double)CGAL::to_double(sqd));
}

float HiMesh::distance_tree(HiMesh *target){
	double min_dist = DBL_MAX;
	for(HiMesh::Vertex_iterator v = vertices_begin(); v != vertices_end(); ++v){
		double dist = target->distance_tree(v->point());
		if(min_dist>dist){
			min_dist = dist;
		}
	}
	for(HiMesh::Vertex_iterator v = target->vertices_begin(); v != target->vertices_end(); ++v){
		double dist = distance_tree(v->point());
		if(min_dist>dist){
			min_dist = dist;
		}
	}
	return sqrt(min_dist);
}

float HiMesh::area(){
	list<Triangle> triangles = get_triangles();
	float a = 0.0;
	for(const Triangle &t:triangles){
		a += tdbase::triangle_area(t);
	}
	return a;
}

float HiMesh::get_volume() {
	Polyhedron *poly = to_triangulated_polyhedron();
	float volume = CGAL::Polygon_mesh_processing::volume(*poly);
	delete poly;
	return volume;
}

/*
 *
 * some utility functions for geometric computations
 *
 * */

float HiMesh::distance(const Point &p1, const Point &p2){
	float dist = 0;
	for(int i=0;i<3;i++){
		dist += (p2[i]-p1[i])*(p2[i]-p1[i]);
	}
	return sqrt(dist);
}

float HiMesh::distance(const Point &p, const Triangle &t){
	return PointTriangleDist((const float *)&p, (const float *)&t);
}

float HiMesh::distance(const tdbase::Point &p, const HiMesh::Face_iterator &fit){
	const float point[3] = {p.x(), p.y(), p.z()};
	HiMesh::Halfedge_const_handle hd = fit->halfedge();
	HiMesh::Halfedge_const_handle h = hd->next();
	float mindist = DBL_MAX;
	while(h->next()!=hd){
		Point p1 = hd->vertex()->point();
		Point p2 = h->vertex()->point();
		Point p3 = h->next()->vertex()->point();
		h = h->next();
		const float triangle[9] = {p1.x(), p1.y(), p1.z(),
								   p2.x(), p2.y(), p2.z(),
								   p3.x(), p3.y(), p3.z()};
		float dist = PointTriangleDist(point, triangle);
		mindist = min(mindist, dist);
	}
	return mindist;
}

float HiMesh::distance(const Triangle &t1, const Triangle &t2){
	return TriDist((const float *)&t1, (const float *)&t2);
}

bool HiMesh::intersect(const Triangle &t1, const Triangle &t2){
	return TriInt((const float *)&t1, (const float *)&t2);
}




}


