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


/**
  * Test if a vertex removal will violate the manifold property of a mesh.
  * \return true if it will else false.
  */
bool MyMesh::willViolateManifold(const std::vector<Halfedge_const_handle> &polygon) const
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


/*
 *
 * Test whether a vertex is protruding
 * Core function for 3DPro
 *
 * */

void MyMesh::profileProtruding(){
	int protruding = 0;
	int recessing = 0;
	for(MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		if(isProtruding(vit)){
			protruding++;
		}else{
			recessing++;
		}
	}
	printf("%d %d %f\n",protruding,recessing,protruding*100.0/(protruding+recessing));
}

bool MyMesh::isProtruding(Vertex_const_handle v) const
{
	  //test convexity
	  std::vector<Halfedge_const_handle> heh_oneRing;
	  heh_oneRing.reserve(v->vertex_degree());
	  //vh_oneRing.push_back(v);
	  Halfedge_around_vertex_const_circulator hit(v->vertex_begin()), end(hit);
	  do
	  {
		  heh_oneRing.push_back(hit->opposite());
	  }
	  while(++hit != end);
	  return isProtruding(heh_oneRing);
}





#define PI 3.1415926
bool MyMesh::isProtruding(const std::vector<Halfedge_const_handle> &polygon) const
{
	if(polygon.size()<2){
		return true;
	}
	Point top(polygon[0]->opposite()->vertex()->point());
	vector<Point> rings;

	for(Halfedge_const_handle h:polygon){
		Point p(h->vertex()->point());
		rings.push_back(p);
	}

	// find one possible triangulation which fulfill protruding
	bool is_recessing = false;
	for(int t=0;t<rings.size()-1;t++){
		is_recessing = false;
		// evaluate all the tetrahedrons
		for(int i=1;i<rings.size()-2;i++){
			// the random vector pointing from the bottom triangle to the top vertex
			double r1 = top.x()-rings[t].x();
			double r2 = top.y()-rings[t].y();
			double r3 = top.z()-rings[t].z();

			// calculate the normal vector of the bottom triangle
			//printf("3 0 %d %d|",i+1,i+2);
			double a1 = rings[t].x()-rings[(i+t)%rings.size()].x();
			double a2 = rings[t].y()-rings[(i+t)%rings.size()].y();
			double a3 = rings[t].z()-rings[(i+t)%rings.size()].z();
			double b1 = rings[t].x()-rings[(i+t+1)%rings.size()].x();
			double b2 = rings[t].y()-rings[(i+t+1)%rings.size()].y();
			double b3 = rings[t].z()-rings[(i+t+1)%rings.size()].z();

			// n*a=0 n*b=0
			double n1 = a2*b3-a3*b2;
			double n2 = a3*b1-a1*b3;
			double n3 = a1*b2-a2*b1;

			// calculate the angle between the normal and vector 1->0
			double cosvalue = (r1*n1+r2*n2+r3*n3)/(sqrt(r1*r1+r2*r2+r3*r3)*sqrt(n1*n1+n2*n2+n3*n3));
			double angle = acos(cosvalue)*180/PI;
			// avoid the border case
			if(angle>90.5){
				is_recessing = true;
			}
			//printf("%d\tangle: %f %f\n",t,cosvalue,angle);
		}
		// this vertex can be protruding
		if(is_recessing == false){
			break;
		}
	}

	// print the removed part into a single polyhedron for visualization
	if(is_recessing && false){
		printf("OFF\n\n");
		printf("%ld %ld 0\n",1+rings.size(),1+rings.size());
		printf("%f %f %f\n",top.x(),top.y(),top.z());
		for(Point p:rings){
			printf("%f %f %f\n",p.x(),p.y(),p.z());
		}
		for(int i=0;i<rings.size()-1;i++){
			printf("3 %d %d 0 0 255 0\n",i+2,i+1);
		}
		printf("3 1 %ld 0 0 255 0\n",rings.size());

		printf("%ld",rings.size());
		for(int i=1;i<=rings.size();i++){
			printf(" %d",i);
		}
		printf(" 255 0 0\n\n");
	}
	// no recessing point
	return !is_recessing;
}

/**
  * Test for the convexity of a polygon
  */
bool MyMesh::isConvex(const std::vector<Vertex_const_handle> & polygon) const
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
  * Test if a polygon is plannar.
  */
bool MyMesh::isPlanar(const std::vector<Vertex_const_handle> &polygon, float epsilon) const
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
float MyMesh::removalError(Vertex_const_handle v,
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


/**
  * Test if a vertex is removable.
  */
bool MyMesh::isRemovable(Vertex_const_handle v) const
{
	if (v != vh_departureConquest[0] && v != vh_departureConquest[1] &&
		!v->isConquered() && v->vertex_degree() > 2 && v->vertex_degree() <= 8)
	{
	  //test convexity
	  std::vector<Vertex_const_handle> vh_oneRing;
	  std::vector<Halfedge_const_handle> heh_oneRing;

	  vh_oneRing.reserve(v->vertex_degree());
	  heh_oneRing.reserve(v->vertex_degree());
	  //vh_oneRing.push_back(v);
	  Halfedge_around_vertex_const_circulator hit(v->vertex_begin()), end(hit);
	  do
	  {
			vh_oneRing.push_back(hit->opposite()->vertex());
			heh_oneRing.push_back(hit->opposite());
	  }
	  while(++hit != end);
	  //&& isConvex(vh_oneRing)
	  bool removable = !willViolateManifold(heh_oneRing)  && isProtruding(heh_oneRing);
	  return removable;
	}
	return false;
}


