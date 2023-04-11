/*
 * himesh_3dpro.cpp
 *
 *  Created on: Apr 10, 2023
 *      Author: teng
 */


#include "himesh.h"

namespace hispeed{

/**
  * Test for the convexity of a polygon
  */
bool HiMesh::isConvex(const std::vector<Vertex_const_handle> & polygon) const
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
bool HiMesh::willViolateManifold(const std::vector<Halfedge_const_handle> &polygon) const
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

bool HiMesh::isProtruding(const std::vector<Halfedge_const_handle> &polygon) const
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

	// find one possible triangulation way that fulfill protruding
	bool is_recessing = false;
	// t is the starting vertex on the ring
	for(int t=0;t<rings.size()-1;t++){
	//for(int t=0;t<1;t++){

		is_recessing = false;
		// evaluate all the tetrahedrons
		for(int i=1;i<rings.size()-1;i++){
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

			if(!global_ctx.counter_clock){
				n1 = -n1;
				n2 = -n2;
				n3 = -n3;
			}

			// calculate the angle between the normal and vector t->0
			double cosvalue = (r1*n1+r2*n2+r3*n3)/(sqrt(r1*r1+r2*r2+r3*r3)*sqrt(n1*n1+n2*n2+n3*n3));
			double angle = acos(cosvalue)*180/M_PI;
			// avoid some corner case, such that bigger than 90.5 instead of 90, increase the tolerance.
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
		printf("%d\n",is_recessing);
		printf("OFF\n");
		printf("%ld %ld 0\n",1+rings.size(),1+rings.size());
		printf("%f %f %f\n",top.x(),top.y(),top.z());
		for(Point p:rings){
			printf("%f %f %f\n",p.x(),p.y(),p.z());
		}

		if(global_ctx.counter_clock){
			for(int i=0;i<rings.size()-1;i++){
				printf("3 0 %d %d 0 255 0\n",i+1,i+2);
			}
			printf("3 0 %ld 1 0 255 0\n",rings.size());

			printf("%ld",rings.size());
			for(int i=rings.size();i>=1;i--){
				printf(" %d",i);
			}
		}else{
			for(int i=0;i<rings.size()-1;i++){
				printf("3 %d %d 0 0 255 0\n",i+2,i+1);
			}
			printf("3 1 %ld 0 0 255 0\n",rings.size());

			printf("%ld",rings.size());
			for(int i=1;i<=rings.size();i++){
				printf(" %d",i);
			}
		}

		printf(" 255 0 0\n");

		vector<Point> norms;

		for(int i=1;i<rings.size()-1;i++){
			// the random vector pointing from the bottom triangle to the top vertex
			double r1 = top.x()-rings[0].x();
			double r2 = top.y()-rings[0].y();
			double r3 = top.z()-rings[0].z();

			// calculate the normal vector of the bottom triangle
			//printf("3 0 %d %d|",i+1,i+2);
			double a1 = rings[0].x()-rings[(i+0)%rings.size()].x();
			double a2 = rings[0].y()-rings[(i+0)%rings.size()].y();
			double a3 = rings[0].z()-rings[(i+0)%rings.size()].z();
			double b1 = rings[0].x()-rings[(i+0+1)%rings.size()].x();
			double b2 = rings[0].y()-rings[(i+0+1)%rings.size()].y();
			double b3 = rings[0].z()-rings[(i+0+1)%rings.size()].z();

			// n*a=0 n*b=0
			double n1 = a2*b3-a3*b2;
			double n2 = a3*b1-a1*b3;
			double n3 = a1*b2-a2*b1;

			if(!global_ctx.counter_clock){
				n1 = -n1;
				n2 = -n2;
				n3 = -n3;
			}
			Point p(n1*100,n2*100,n3*100);
			norms.push_back(p);
			//printf("%d\tangle: %f %f\n",t,cosvalue,angle);
		}

		// the vertical lines
		printf("OFF\n");
		printf("%ld %ld 0\n",2*(rings.size()-2)+1,rings.size()-2);
		for(int i=1;i<rings.size()-1;i++){
			printf("%f %f %f\n",rings[i].x(),rings[i].y(),rings[i].z());
		}
		for(int i=1;i<rings.size()-1;i++){
			printf("%f %f %f\n",rings[i].x()+norms[i-1].x(),rings[i].y()+norms[i-1].y(),rings[i].z()+norms[i-1].z());
		}
		printf("%f %f %f\n",rings[0].x(),rings[0].y(),rings[0].z());
		for(int i=1;i<rings.size()-1;i++){
			if(i==1){
				printf("3 %d %ld %ld 255 0 0\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}else if(i==2){
				printf("3 %d %ld %ld 0 255 0\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}else if(i==3){
				printf("3 %d %ld %ld 0 0 255\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}else{
				printf("3 %d %ld %ld\n",i-1,i-1+rings.size()-2,2*(rings.size()-2));
			}
		}
		exit(0);
	}
	// no recessing point
	return !is_recessing;
}

/**
  * Test if a vertex is removable.
  */
bool HiMesh::isRemovable(Vertex_handle v) const
{
//	if(size_of_vertices()<10){
//		return false;
//	}
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
	  //
	  bool removable = !willViolateManifold(heh_oneRing);
	  //&& isConvex(vh_oneRing);
	  if(removable && !isProtruding(heh_oneRing)){
		 v->setRecessing();
	  }
	  return removable;
	}
	return false;
}

/*
 *
 * Test whether a vertex is protruding
 * Core function for 3DPro
 *
 * */

void HiMesh::profileProtruding(){
	int protruding = 0;
	int recessing = 0;
	for(HiMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit){
		std::vector<Halfedge_const_handle> heh_oneRing;
		heh_oneRing.reserve(vit->vertex_degree());
		//vh_oneRing.push_back(v);
		Halfedge_around_vertex_const_circulator hit(vit->vertex_begin()), end(hit);
		do
		{
		  heh_oneRing.push_back(hit->opposite());
		}while(++hit != end);
		if(isProtruding(heh_oneRing)){
			protruding++;
		}else{
			recessing++;
		}
	}
	printf("%d %d %f\n",protruding,recessing,protruding*100.0/(protruding+recessing));
}

pair<float, float> HiMesh::getHausdorfDistance(){
	assert(i_nbDecimations>=i_curDecimationId);
	return i_nbDecimations>i_curDecimationId?globalHausdorfDistance[i_nbDecimations - i_curDecimationId-1]:std::pair<float, float>(0, 0);
}

pair<float, float> HiMesh::getNextHausdorfDistance(){
	assert(i_nbDecimations>i_curDecimationId);
	return i_nbDecimations>(i_curDecimationId+1)?globalHausdorfDistance[i_nbDecimations - i_curDecimationId - 2]:std::pair<float, float>(0, 0);
}

static float point_to_face_distance(const Point &p, const HiMesh::Face_iterator &fit){
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

void HiMesh::computeHausdorfDistance(){

	pair<float, float> current_hausdorf = pair<float, float>(0.0, 0.0);
	float dist = DBL_MAX;

	struct timeval start = get_cur_time();
	get_aabb_tree_triangle();
	double hdist = 0;
	for(Point &p:removedPoints){
		double dist = this->distance_tree(p);
		//triangle_tree->closest_point_and_primitive(query, hint)
		hdist = max(hdist, dist);
	}
	clear_aabb_tree();
	//logt("%f %d", start, hdist, removedPoints.size());

	log("%ld groups", this->map_group.size());

	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
		if(fit->isSplittable()){
			vector<Point> &ps = fit->getImpactPoints();
			if(ps.size()==0){
			  continue;
			}
			float maxdist = 0.0;
			for(Point &p:ps){
				float dist = point_to_face_distance(p, fit);
				maxdist = max(dist, maxdist);
			}
			//log("%f", farthest_point);
			fit->setProgressive(maxdist);
			fit->setConservative(0.0);
			//log("%f", fit->getHausdorfDistance().second);
		}
		current_hausdorf.second = max(current_hausdorf.second, fit->getHausdorfDistance().second);
		dist = min(dist, fit->getHausdorfDistance().second);
	}
	//logt("%f", start, current_hausdorf.second);
	globalHausdorfDistance.push_back(current_hausdorf);
	if(global_ctx.verbose>=2){
		log("encode %d:\t[%.2f %.2f]\t%ld", i_curDecimationId, dist, current_hausdorf.second, size_of_vertices());
	}
}


}
