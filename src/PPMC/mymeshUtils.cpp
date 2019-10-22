/*****************************************************************************
* Copyright (C) 2011 Adrien Maglo
*
* This file is part of PPMC.
*
* PPMC is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* PPMC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with PPMC.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include "../PPMC/mymesh.h"


Vector MyMesh::computeNormal(Facet_const_handle f) const
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


Vector MyMesh::computeNormal(Halfedge_const_handle heh_gate) const
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


Vector MyMesh::computeNormal(const std::vector<Vertex_const_handle> & polygon) const
{
    Vector n(0,0,0);
    int s = polygon.size();
    for(int i=0; i<s; ++i)
    {
        Vector op(CGAL::ORIGIN, polygon[i]->point());
        
     //    std::cerr << "polygon " << i << "\t" <<  polygon[i]->point()  << std::endl;
        Vector op2(CGAL::ORIGIN, polygon[(i+1)%s]->point());
        n = n + CGAL::cross_product(op,op2);
    }
    float f_sqLen = n.squared_length();
 //   std::cerr << "f_sqLen" << f_sqLen << std::endl;
 //   exit(0);
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


// Compute the normal of a triangle.
Vector MyMesh::computeNormal(Point p[3]) const
{
    Vector n = CGAL::cross_product(p[1] - p[0], p[2] - p[1]);
    float f_sqLen = n.squared_length();
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


Vector MyMesh::computeVertexNormal(Halfedge_const_handle heh) const
{
    MyMesh::Halfedge_around_vertex_const_circulator hit = heh->vertex_begin(), hit_end = hit;
    Vector n(CGAL::NULL_VECTOR);
    CGAL_For_all(hit, hit_end)
    {
        n = n + computeNormal(hit->opposite());
    }
    float f_sqLen = n.squared_length();
    return f_sqLen == 0 ? CGAL::NULL_VECTOR : n / sqrt(f_sqLen);
}


Point MyMesh::barycenter(Facet_const_handle f) const
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


Point MyMesh::barycenter(Halfedge_handle heh_gate) const
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


Point MyMesh::barycenter(const std::vector<Vertex_const_handle> &polygon) const
{
    Point b = CGAL::ORIGIN;
    unsigned i_size = polygon.size();
    for (unsigned i = 0; i < i_size; ++i)
          b = b + Vector(CGAL::ORIGIN, polygon[i]->point()) / i_size;
    return b;
}


/**
  * Return the degree of a vertex without counting the new edges.
  */
unsigned MyMesh::vertexDegreeNotNew(Vertex_const_handle vh) const
{
    unsigned i_degree = 0;

    Halfedge_around_vertex_const_circulator hit = vh->vertex_begin(), hit_end = hit;
    CGAL_For_all(hit, hit_end)
    {
        if (!hit->isNew())
            i_degree++;
    }

    return i_degree;
}


/**
  * Return the average laplacian vector of the vertices of a face.
  */
VectorInt MyMesh::avgLaplacianVect(Halfedge_handle heh_gate) const
{
    Halfedge_handle heh = heh_gate;
    VectorInt sumLaplacian = CGAL::NULL_VECTOR;

    do
    {
        VectorInt laplacian = CGAL::NULL_VECTOR;

        Halfedge_handle heh_next = heh->next();
        Halfedge_handle heh2 = heh_next;
        do
        {
            assert(!heh2->isNew() && !heh->isNew());
            laplacian = laplacian + (getQuantizedPos(heh2->vertex()->point())
                                     - getQuantizedPos(heh->vertex()->point()));

            // Move to next neighbor vertex not new.
            do
                heh2 = heh2->prev()->opposite();
            while (heh2->isNew());
        }
        while(heh2 != heh_next);

        sumLaplacian = sumLaplacian + laplacian / vertexDegreeNotNew(heh->vertex());

        // Move to next vertex of the face.
        heh = heh->next();
    }
    while (heh != heh_gate);

    return sumLaplacian / heh_gate->facet_degree();
}


/**
  * Compute the surface of a triangle using Heron's formula.
  */
float MyMesh::triangleSurface(const Point p[]) const
{
    float a = sqrt((p[0] - p[1]).squared_length());
    float b = sqrt((p[1] - p[2]).squared_length());
    float c = sqrt((p[2] - p[0]).squared_length());
    float s = (a + b + c) / 2;

    float f_sqSurface = s * (s - a) * (s - b) * (s - c);

    return f_sqSurface <= 0 ? 0 : sqrt(f_sqSurface);
}

/**
  * Returns an edge length.
  */
float MyMesh::edgeLen(Halfedge_const_handle heh) const
{
    return sqrt((heh->vertex()->point()
                 - heh->opposite()->vertex()->point()).squared_length());
}


/**
  * Returns a face perimter.
  */
float MyMesh::facePerimeter(const Face_handle fh) const
{
    float f_ret = 0;
    Halfedge_around_facet_const_circulator hit = fh->facet_begin(), hit_end = hit;
    CGAL_For_all(hit, hit_end)
    {
        f_ret += edgeLen(hit);
    }
    return f_ret;
}


/**
  * Returns a face surface.
  */
float MyMesh::faceSurface(Halfedge_handle heh) const
{
    unsigned i_degree = heh->facet_degree();
    Point *polygonPos = new Point[i_degree];

    // Store all the face vertices.
    Halfedge_handle hIt = heh;
    unsigned i = 0;
    do
    {
        polygonPos[i++] = hIt->vertex()->point();
        hIt = hIt->next();
    }
    while (hIt != heh);

    Point b = barycenter(heh);

    // Compute the surface.
    float f_ret = 0;
    for (i = 0; i < i_degree; i++)
    {
        Point trianglePos[3];
        trianglePos[0] = b;
        trianglePos[1] = polygonPos[i];
        trianglePos[2] = polygonPos[(i + 1) % i_degree];
        f_ret += triangleSurface(trianglePos);
    }

    assert(f_ret >= 0);

    delete[] polygonPos;
    return f_ret;
}


/**
  * Push the first halfedge for the coding and decoding conquest in the gate queue.
  */
void MyMesh::pushHehInit()
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
