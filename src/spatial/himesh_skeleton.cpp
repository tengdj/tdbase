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
				double dist = distance(vertices[n], ret[c]);
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
 * around the points of the skeleton, and generate a list
 * of axis aligned boxes for those sets
 * */
vector<Voxel *> HiMesh::generate_voxels_skeleton(int voxel_size){
	timeval start = tdbase::get_cur_time();

	int voxel_num = size_of_facets()/voxel_size+1;
	//cout << voxel_num <<" "<< size_of_facets() << " " << voxel_size << endl;
	
	vector<Voxel *> voxels;
	aab box = get_mbb();
	// avoid using too few voxels
	if(voxel_num<=3){
		Voxel *v = new Voxel();
		v->set_box(box);
		voxels.push_back(v);
		return voxels;
	}

	vector<Point> skeleton_points = get_skeleton_points(voxel_num);
	for(int i=0;i<skeleton_points.size();i++){
		Voxel *v = new Voxel();
		v->core[0] = skeleton_points[i][0];
		v->core[1] = skeleton_points[i][1];
		v->core[2] = skeleton_points[i][2];
		voxels.push_back(v);
	}

	// assigning each facet to proper voxel
	for (Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		Point p = HiMesh::barycenter(f);
		aab box = HiMesh::bounding_box(f);
		float min_dist = DBL_MAX;
		int gid = -1;
		for(int j=0;j<skeleton_points.size();j++){
			float cur_dist = distance(p, skeleton_points[j]);
			if(cur_dist<min_dist){
				gid = j;
				min_dist = cur_dist;
			}
		}
		voxels[gid]->update(box);
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

/*
 *
 * conduct voxelization
 *
 * */
vector<Voxel *> HiMesh::voxelization(int voxel_size){
	vector<Voxel *> voxels;
	if(voxel_size<=1){
		Voxel *vox = new Voxel();
		vox->set_box(get_mbb());
		voxels.push_back(vox);
	}

	aab box = get_mbb();
	float min_dim = std::min(box.high[2]-box.low[2], std::min(box.high[1]-box.low[1], box.high[0]-box.low[0]));
	float div = (box.high[2]-box.low[2])*(box.high[1]-box.low[1])*(box.high[0]-box.low[0])/(min_dim*min_dim*min_dim);
	float multi = std::pow(1.0*voxel_size/div, 1.0/3);

	int dim[3];
	for(int i=0;i<3;i++){
		dim[i] = ((box.high[i]-box.low[i])*multi/min_dim+0.5);
		assert(dim[i]>0);
	}

	log("%d %d %d",dim[0],dim[1],dim[2]);

	bool *taken = new bool[dim[0]*dim[1]*dim[2]];
	for(int i=0;i<dim[0]*dim[1]*dim[2];i++){
		taken[i] = false;
	}

	unordered_set<Point> points;
	int old_sampled_rate = sampling_rate;
	sampling_rate = 40;
	sample_points(points);
	sampling_rate = old_sampled_rate;

	for(Point p:points){
		int x = (p.x()-box.low[0])*dim[0]/(box.high[0]-box.low[0]);
		int y = (p.y()-box.low[1])*dim[1]/(box.high[1]-box.low[1]);
		int z = (p.z()-box.low[2])*dim[2]/(box.high[2]-box.low[2]);

		if(x==dim[0]){
			x = dim[0]-1;
		}
		if(y==dim[1]){
			y = dim[1]-1;
		}
		if(z==dim[2]){
			z = dim[2]-1;
		}
		assert(x<dim[0] && y<dim[1] && z<dim[2]);

		int idx = z*dim[1]*dim[0]+y*dim[0]+x;
		if(!taken[idx]){
			Voxel *vox = new Voxel();
			vox->low[0] = x*(box.high[0]-box.low[0])/dim[0];
			vox->low[1] = y*(box.high[1]-box.low[1])/dim[1];
			vox->low[2] = z*(box.high[2]-box.low[2])/dim[2];
			vox->high[0] = (x+1)*(box.high[0]-box.low[0])/dim[0];
			vox->high[1] = (y+1)*(box.high[1]-box.low[1])/dim[1];
			vox->high[2] = (z+1)*(box.high[2]-box.low[2])/dim[2];
			voxels.push_back(vox);
		}
		taken[idx] = true;
	}
	return voxels;
}

}
