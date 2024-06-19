/*
 * himesh_utils.cpp
 *
 *  Created on: Apr 6, 2023
 *      Author: teng
 */

#include "himesh.h"

namespace tdbase{

Vector HiMesh::computeNormal(Facet_const_handle f)
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


Vector HiMesh::computeNormal(Halfedge_const_handle heh_gate)
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


Vector HiMesh::computeNormal(const std::vector<Vertex_const_handle> & polygon)
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
Vector HiMesh::computeNormal(Point p[3])
{
    Vector n = CGAL::cross_product(p[1] - p[0], p[2] - p[1]);
    float f_sqLen = n.squared_length();
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


Vector HiMesh::computeVertexNormal(Halfedge_const_handle heh)
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


Point HiMesh::barycenter(Facet_const_handle f) {
  //compute the barycenter of the face:
  Point barycenter = CGAL::ORIGIN;
  double d = f->facet_degree();
  Halfedge_around_facet_const_circulator hit(f->facet_begin()), end(hit);
  do {
	  barycenter = barycenter + Vector(CGAL::ORIGIN, hit->vertex()->point())/d;
  }
  while(++hit != end);
  return barycenter;
}

aab HiMesh::bounding_box(Facet_const_handle f){
	aab box;
	Halfedge_around_facet_const_circulator hit(f->facet_begin()), end(hit);
	do {
		Point p = hit->vertex()->point();
		box.update(p.x(), p.y(), p.z());
	}while(++hit != end);
	return box;
}

Point HiMesh::barycenter(Halfedge_handle heh_gate) {
    Halfedge_handle heh = heh_gate;
    Vector barycenter = CGAL::NULL_VECTOR;
    unsigned i_degree = 0;

    do {
        barycenter = barycenter + Vector(CGAL::ORIGIN, heh->vertex()->point());
        heh = heh->next();
        i_degree++;
    }
    while (heh != heh_gate);

    return CGAL::ORIGIN + barycenter / i_degree;
}


Point HiMesh::barycenter(const std::vector<Vertex_const_handle> &polygon) {
    Point b = CGAL::ORIGIN;
    unsigned i_size = polygon.size();
    for (unsigned i = 0; i < i_size; ++i)
          b = b + Vector(CGAL::ORIGIN, polygon[i]->point()) / i_size;
    return b;
}

/**
  * Compute the surface of a triangle using Heron's formula.
  */
float HiMesh::triangleSurface(const Point p[]) {
    float a = sqrt((p[0] - p[1]).squared_length());
    float b = sqrt((p[1] - p[2]).squared_length());
    float c = sqrt((p[2] - p[0]).squared_length());
    float s = (a + b + c) / 2;

    float f_sqSurface = s * (s - a) * (s - b) * (s - c);

    return f_sqSurface <= 0 ? 0 : sqrt(f_sqSurface);
}

/**
  * Test if a polygon is plannar.
  */
bool HiMesh::isPlanar(const std::vector<Vertex_const_handle> &polygon, float epsilon)
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
  * Test for the convexity of a polygon
  */
bool HiMesh::isConvex(const std::vector<Vertex_const_handle> & polygon)
{
  //project all points on a plane, taking the first point as origin
  Vector n = computeNormal(polygon);
  int s = polygon.size();
  std::vector<Point> projPoints(s);
  for(int i=0; i<s; ++i)
  {
        //project polygon[i]->point() on the plane with normal n
        projPoints[i] = polygon[i]->point() - n*(Vector(polygon[0]->point(), polygon[i]->point())*n);
  }

  //now use the following test: a polygon is concave if for each edge, all the other points lie on the same side of the edge
  for(int i=0; i<s; ++i)
  {
        Vector ev(projPoints[i], projPoints[(i+1)%s]);
        int globalSide = 0;
        int comp[9] = {0,1,2,1,1,3,2,3,2};
        //(0,0) -> 0
        //(0,+) -> +
        //(0,-) -> -
        //(+,0) -> +
        //(+,+) -> +
        //(+,-) -> 3
        //(-,0) -> -
        //(-,+) -> 3
        //(-,-) -> -
        for(int j=0; j<s; ++j)
        {
          if( j==i || j==(i+1) )
                continue;
          Vector dv(projPoints[i], projPoints[j]);
          Vector evxn = CGAL::cross_product(ev,n);
          double cp = evxn*dv;
          int side = (fabs(cp) > 0.000001)*(cp>0 ? 1 : 2);
          globalSide = comp[globalSide*3+side];
          if(globalSide==3)
          {
        	  	// non convex
                return false;
          }
        }
  }
  //convex
  return true;
}

/**
  * Test if a vertex removal will violate the manifold property of a mesh.
  * \return true if it will else false.
  */
bool HiMesh::willViolateManifold(const std::vector<Halfedge_const_handle> &polygon)
{
    unsigned i_degree = polygon.size();

    // Test if a patch vertex is not connected to one vertex
    // that is not one of its direct neighbor.
    // Test also that two vertices of the patch will not be doubly connected
    // after the vertex cut opeation.
    for (unsigned i = 0; i < i_degree; ++i)
    {
        Halfedge_around_vertex_const_circulator Hvc = polygon[i]->vertex()->vertex_begin();
        Halfedge_around_vertex_const_circulator Hvc_end = Hvc;
        CGAL_For_all(Hvc, Hvc_end)
        {
            // Look if the current vertex belongs to the patch.
            Vertex_const_handle vh = Hvc->opposite()->vertex();
            for (unsigned j = 0; j < i_degree; ++j)
            {
                if (vh == polygon[j]->vertex())
                {
                    unsigned i_prev = i == 0 ? i_degree - 1 : i - 1;
                    unsigned i_next = i == i_degree - 1 ? 0 : i + 1;

                    if ((j == i_prev && polygon[i]->facet_degree() != 3) // The vertex cut operation is forbidden.
                        || (j == i_next && polygon[i]->opposite()->facet_degree() != 3)) // The vertex cut operation is forbidden.
                        return true;
                }
            }
        }
    }

    return false;
}

/**
  * Compute the error introduced by the removal of a vertex.
  * Same metric as in the progressive coder of Alliez and Desbrun.
  */
float HiMesh::removalError(Vertex_const_handle v, const std::vector<Vertex_const_handle> &polygon) {
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

}

