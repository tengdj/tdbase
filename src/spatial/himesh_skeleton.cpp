/*
 * himesh_internal_index.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */

#include "himesh.h"

namespace tdbase{

/*
 *
 * extract the skeleton of current polyhedron
 *
 * */
vector<Point> HiMesh::get_skeleton_points(int num_cores){
	vector<Point> ret;
	Point *vertices = NULL;

	const size_t numv = fill_vertices((float *&)vertices);
	const size_t gap = numv/num_cores + 1;
	for(size_t n=0;n<numv;n+=gap){
		ret.push_back(vertices[n]);
	}

	Point *tmp = new Point[ret.size()];
	uint *ct = new uint[ret.size()];
	Point zerop(0.0, 0.0, 0.0);

	for(int i=0;i<ret.size();i++){
		tmp[i] = zerop;
		ct[i] = 0;
	}

	struct timeval start = get_cur_time();
	for(int it=0;it<2;it++){
		for(size_t n=0;n<numv;n++){
			size_t closest = 0;
			double cur_min = DBL_MAX;
			for(size_t c=0;c<ret.size();c++){
				double dist = get_distance(vertices[n], ret[c]);
				if(dist < cur_min){
					closest = c;
					cur_min = dist;
				}
			}
			for(int d=0;d<3;d++){
				((float *)&tmp[closest])[d] += vertices[n][d];
			}
			ct[closest]++;
		}

		for(size_t c=0;c<ret.size();c++){
			if(ct[c]>0){
				for(int d=0;d<3;d++){
					((float *)&ret[c])[d] = tmp[c][d]/ct[c];
				}
			}
		}
	}

	delete []vertices;
	return ret;
}

/*
 * get the skeleton of the polyhedron, and then group the triangles
 * or edges around the points of the skeleton, and generate a list
 * of axis aligned boxes for those sets
 * */
vector<Voxel *> HiMesh::generate_voxels_skeleton(int voxel_num){

	voxel_num = std::max(1, voxel_num);
	timeval start = tdbase::get_cur_time();
	vector<Voxel *> voxels;
	aab box = get_mbb();
	if(voxel_num<=1){
		Voxel *v = new Voxel();
		v->set_box(box);
		voxels.push_back(v);
		return voxels;
	}

	//log("%d %d ",size_of_vertices(), voxel_num);
	//assert(num_cores>3);
	// this step takes 99 percent of the computation load
	vector<Point> skeleton_points = get_skeleton_points(voxel_num);
	for(int i=0;i<skeleton_points.size();i++){
		Voxel *v = new Voxel();
		v->core[0] = skeleton_points[i][0];
		v->core[1] = skeleton_points[i][1];
		v->core[2] = skeleton_points[i][2];
		voxels.push_back(v);
	}

	// return one single box if less than 2 points are sampled
	if(voxels.size()==0){
		Voxel *v = new Voxel();
		v->set_box(box);
		voxels.push_back(v);
		return voxels;
	}else if(voxels.size()==1){
		voxels[0]->set_box(box);
		return voxels;
	}

	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		Point p1 = f->halfedge()->vertex()->point();
		Point p2 = f->halfedge()->next()->vertex()->point();
		Point p3 = f->halfedge()->next()->next()->vertex()->point();
		Point p((p1[0]+p2[0]+p3[0])/3, (p1[1]+p2[1]+p3[1])/3,(p1[2]+p2[2]+p3[2])/3);
		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<skeleton_points.size();j++){
			float cur_dist = 0;
			for(int i=0;i<3;i++){
				cur_dist += (skeleton_points[j][i]-p[i])*(skeleton_points[j][i]-p[i]);
			}
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		voxels[gid]->update(p1.x(),p1.y(),p1.z());
		voxels[gid]->update(p2.x(),p2.y(),p2.z());
		voxels[gid]->update(p3.x(),p3.y(),p3.z());

		// we do not need the facets in each voxel, such that we do not need to insert one facet
		voxels[gid]->num_triangles++;
	}

	// erase the one without any data in it
	int vs=voxels.size();
	for(int i=0;i<vs;){
		if(voxels[i]->num_triangles==0){
			voxels.erase(voxels.begin()+i);
			vs--;
		}else{
			// reset to 0 as we never truly inserted any facets
			voxels[i]->num_triangles = 0;
			i++;
		}
	}
	skeleton_points.clear();

	return voxels;
}

}
