/*
 * himesh_utils.cpp
 *
 *  Created on: Apr 6, 2023
 *      Author: teng
 */

#include "himesh.h"

namespace hispeed{

Vector HiMesh::computeNormal(Facet_const_handle f) const
{
    Halfedge_around_facet_const_circulator hit(f->facet_begin()), end(hit);
    Vector n(0,0,0);
    do
    {
        Vector op(CGAL::ORIGIN, hit->vertex()->point());
        Vector op2(CGAL::ORIGIN, hit->next()->vertex()->point());
        n = n + CGAL::cross_product(op,op2);
    }
    while(++hit != end);

    float f_sqLen = n.squared_length();
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


Vector HiMesh::computeNormal(Halfedge_const_handle heh_gate) const
{
    Halfedge_const_handle heh = heh_gate;
    std::vector<Vertex_const_handle> polygon;
    polygon.reserve(10);

    do
    {
        polygon.push_back(heh->vertex());
        heh = heh->next();
    }
    while (heh != heh_gate);

    return computeNormal(polygon);
}


Vector HiMesh::computeNormal(const std::vector<Vertex_const_handle> & polygon) const
{
    Vector n(0,0,0);
    int s = polygon.size();
    for(int i=0; i<s; ++i)
    {
        Vector op(CGAL::ORIGIN, polygon[i]->point());
        Vector op2(CGAL::ORIGIN, polygon[(i+1)%s]->point());
        n = n + CGAL::cross_product(op,op2);
    }
    float f_sqLen = n.squared_length();
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


// Compute the normal of a triangle.
Vector HiMesh::computeNormal(Point p[3]) const
{
    Vector n = CGAL::cross_product(p[1] - p[0], p[2] - p[1]);
    float f_sqLen = n.squared_length();
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


Vector HiMesh::computeVertexNormal(Halfedge_const_handle heh) const
{
    HiMesh::Halfedge_around_vertex_const_circulator hit = heh->vertex_begin(), hit_end = hit;
    Vector n(CGAL::NULL_VECTOR);
    CGAL_For_all(hit, hit_end)
    {
        n = n + computeNormal(hit->opposite());
    }
    float f_sqLen = n.squared_length();
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


Point HiMesh::barycenter(Facet_const_handle f) const
{
   //compute the barycenter of the face:
  Point barycenter = CGAL::ORIGIN;
  double d = f->facet_degree();
  Halfedge_around_facet_const_circulator hit(f->facet_begin()), end(hit);
  do
  {
        barycenter = barycenter + Vector(CGAL::ORIGIN, hit->vertex()->point())/d;
  }
  while(++hit != end);
  return barycenter;
}


Point HiMesh::barycenter(Halfedge_handle heh_gate) const
{
    Halfedge_handle heh = heh_gate;
    Vector barycenter = CGAL::NULL_VECTOR;
    unsigned i_degree = 0;

    do
    {
        barycenter = barycenter + Vector(CGAL::ORIGIN, heh->vertex()->point());
        heh = heh->next();
        i_degree++;
    }
    while (heh != heh_gate);

    return CGAL::ORIGIN + barycenter / i_degree;
}


Point HiMesh::barycenter(const std::vector<Vertex_const_handle> &polygon) const
{
    Point b = CGAL::ORIGIN;
    unsigned i_size = polygon.size();
    for (unsigned i = 0; i < i_size; ++i)
          b = b + Vector(CGAL::ORIGIN, polygon[i]->point()) / i_size;
    return b;
}

/**
  * Compute the surface of a triangle using Heron's formula.
  */
float HiMesh::triangleSurface(const Point p[]) const
{
    float a = sqrt((p[0] - p[1]).squared_length());
    float b = sqrt((p[1] - p[2]).squared_length());
    float c = sqrt((p[2] - p[0]).squared_length());
    float s = (a + b + c) / 2;

    float f_sqSurface = s * (s - a) * (s - b) * (s - c);

    return f_sqSurface <= 0 ? 0 : sqrt(f_sqSurface);
}


/**
  * Push the first halfedge for the coding and decoding conquest in the gate queue.
  */
void HiMesh::pushHehInit()
{
    // Find the first halfedge.
    Halfedge_handle hehBegin;
    Halfedge_around_vertex_circulator hit(vh_departureConquest[0]->vertex_begin());
    while (1)
    {
        hehBegin = hit->opposite();
        if (hehBegin->vertex() == vh_departureConquest[1])
            break;
        ++hit;
    }
    // Push it to the queue.
    gateQueue.push(hehBegin);
}

/**
  * Test if a polygon is plannar.
  */
bool HiMesh::isPlanar(const std::vector<Vertex_const_handle> &polygon, float epsilon) const
{
    Point b = barycenter(polygon);
    unsigned s = polygon.size();

    //compute polygon caracteristic size
    float caracteristicSize = 0;
    for (unsigned i = 0; i < s; ++i)
        caracteristicSize += Vector(polygon[i]->point(), b).squared_length();
    caracteristicSize/=s;

    Vector n = computeNormal(polygon);
    for (unsigned i = 0; i < s; ++i)
    {
        float dist = n*Vector(b, polygon[i]->point());
        if (dist*dist > epsilon*caracteristicSize)
            return false;
    }
    return true;
}


/**
  * Compute the error introduced by the removal of a vertex.
  * Same metric as in the progressive coder of Alliez and Desbrun.
  */
float HiMesh::removalError(Vertex_const_handle v,
                           const std::vector<Vertex_const_handle> &polygon) const
{
    // Compute the polygon barycenter.
    Point b = barycenter(polygon);

    // Compute the polygon perimeter.
    float f_perimeter = 0;
    unsigned i_size = polygon.size();
    for (unsigned i = 0; i < i_size; ++i)
        f_perimeter += sqrt((polygon[(i + 1) % i_size]->point() - polygon[i]->point()).squared_length());

    // Compute the polygon area.
    float f_area = 0;
    for (unsigned i = 0; i < i_size; i++)
    {
        Point p[3];
        p[0] = polygon[i]->point();
        p[1] = polygon[(i + 1) % i_size]->point();
        p[2] = b;
        f_area += triangleSurface(p);
    }

    // Compute the volume of the polyhedron by using the divergence theorem.
    // http://en.wikipedia.org/wiki/Polyhedron#Volume
    // The polyhedron is in our case a kind of pyramid with the polygon as basis.
    float f_volume = 1 / 3.0 * f_area * Vector(CGAL::ORIGIN, b) * -computeNormal(polygon); // Polygon face.

    // Triangle faces.
    for (unsigned i = 0; i < i_size; i++)
    {
        Point p[3];
        p[0] = polygon[i]->point();
        p[1] = polygon[(i + 1) % i_size]->point();
        p[2] = b;
        Vector n = computeNormal(p);
        f_volume += 1 / 3.0 * triangleSurface(p) * Vector(CGAL::ORIGIN, p[0]) * n;
    }

    float f_metric = f_volume > 0 ? pow(f_volume, 1.0 / 3) * i_size / f_perimeter : 0;

    return f_metric;
}

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

