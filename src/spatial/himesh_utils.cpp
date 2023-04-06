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

}

