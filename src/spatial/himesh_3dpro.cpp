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
	  // && isProtruding(heh_oneRing);
	  //&& isConvex(vh_oneRing)
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

float HiMesh::getHausdorffDistance(){
	assert(i_nbDecimations>=i_curDecimationId);
	return i_nbDecimations>i_curDecimationId?globalHausdorfDistance[i_nbDecimations - i_curDecimationId-1].second:0;
}

float HiMesh::getProxyHausdorffDistance(){
	assert(i_nbDecimations>=i_curDecimationId);
	return i_nbDecimations>i_curDecimationId?globalHausdorfDistance[i_nbDecimations - i_curDecimationId-1].first:0;
}

//pair<float, float> HiMesh::getNextHausdorfDistance(){
//	assert(i_nbDecimations>i_curDecimationId);
//	return i_nbDecimations>(i_curDecimationId+1)?globalHausdorfDistance[i_nbDecimations - i_curDecimationId - 2]:std::pair<float, float>(0, 0);
//}

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

uint HiMesh::sampling_rate = 3;
int HiMesh::calculate_method = 3;

void sample_points_triangle(const Triangle &tri, unordered_set<Point> &points, int num_points){
	const Point &p1 = tri[0];
	const Point &p2 = tri[1];
	const Point &p3 = tri[2];
	points.emplace(p1);
	points.emplace(p2);
	points.emplace(p3);
	if(num_points>3){
		assert(num_points>3);
		int dimx = sqrt(num_points-3);
		int dimy = dimx==0?0:(num_points-3+dimx-1)/dimx;

		Point v1 = p1;
		Point v2(p2.x()-p1.x(), p2.y()-p1.y(), p2.z()-p1.z());
		Point v3(p3.x()-p1.x(), p3.y()-p1.y(), p3.z()-p1.z());

		float step_x = 1.0/(dimx+1);
		float step_y = 1.0/(dimy+1);
		for(float u = 0;u<1;u += step_x){
			for(float v = 0;v<1-u;v += step_y){
				if(!((u==0&&v==0)||(u==1&&v==1))){
					points.emplace(Point(v1.x()+u*v2.x()+v*v3.x(), v1.y()+u*v2.y()+v*v3.y(), v1.z()+u*v2.z()+v*v3.z()));
				}
			}
		}
	}
}

void HiMesh::sample_points(const Triangle &tri, unordered_set<Point> &points, float area_unit){
	int num_points = triangle_area(tri)/area_unit;
	sample_points_triangle(tri, points, num_points);
}

Triangle max_tri;
Triangle min_tri;
float max_area = 0;
float min_area = DBL_MAX;

void HiMesh::sample_points(const HiMesh::Face_iterator &fit, unordered_set<Point> &points, float area_unit){

	const auto hd = fit->halfedge();
	auto h = hd->next();
	while(h->next()!=hd){
		Point p1 = hd->vertex()->point();
		Point p2 = h->vertex()->point();
		Point p3 = h->next()->vertex()->point();
		Triangle tri(p1, p2, p3);
		int num_points = triangle_area(tri)/area_unit+1;
		//cout<<num_points<<endl;
		//num_points = 10;
		int ori = points.size();
		sample_points_triangle(tri, points, num_points);
//		if(num_points < 3){
//			cout<<triangle_area(tri)<<" "<<area_unit<<" "<<num_points<<" "<<points.size()-ori<<endl;
//		}
//		if(max_area<triangle_area(tri)){
//			max_area = triangle_area(tri);
//			max_tri = tri;
//		}
//		if(min_area > triangle_area(tri)){
//			min_area = triangle_area(tri);
//			min_tri = tri;
//		}
		h = h->next();
	}
}

void HiMesh::sample_points(unordered_set<Point> &points){
	const float area_unit = get_sample_density();

	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
		sample_points(fit, points, area_unit);
	}
	log("%.5f %.10f (%d*%d)=%d %d", area(), area_unit, size_of_triangles(), sampling_rate, size_of_triangles()*sampling_rate, points.size());

//	points.clear();
//	sample_points_triangle(max_tri, points, triangle_area(max_tri)/area_unit);
//	sample_points_triangle(min_tri, points, 10);
//	cout<<min_tri<<endl;
//	log("%d %d", (int)(triangle_area(max_tri)/(step*step)), (int)(a/triangle_area(min_tri)));
//	log("%.15f %.15f %.15f %.15f %d", triangle_area(max_tri)/(step*step), triangle_area(min_tri), area_unit, a, total_num_points);
}

Triangle expand(Triangle &tri){
	const Point &p1 = tri[0];
	const Point &p2 = tri[1];
	const Point &p3 = tri[2];
	return Triangle(Point(2*p1.x()-p2.x()/2-p3.x()/2, 2*p1.y()-p2.y()/2-p3.y()/2, 2*p1.z()-p2.z()/2-p3.z()/2),
					Point(2*p2.x()-p1.x()/2-p3.x()/2, 2*p2.y()-p1.y()/2-p3.y()/2, 2*p2.z()-p1.z()/2-p3.z()/2),
					Point(2*p3.x()-p1.x()/2-p2.x()/2, 2*p3.y()-p1.y()/2-p2.y()/2, 2*p3.z()-p1.z()/2-p2.z()/2));
}

vector<Triangle> triangulate(HiMesh::Face_iterator fit){
	vector<Triangle> ret;
	const auto hd = fit->halfedge();
	auto h = hd->next();
	while(h->next()!=hd){
		Point p1 = hd->vertex()->point();
		Point p2 = h->vertex()->point();
		Point p3 = h->next()->vertex()->point();
		ret.push_back(Triangle(p1,p2,p3));
		h = h->next();
	}
	return ret;
}

int tt = 0;

// TODO: a critical function, need to be further optimized
void HiMesh::computeHausdorfDistance(){
	// do not compute
	if(calculate_method == HCT_NULL ){
		globalHausdorfDistance.push_back({0, 0});
		return;
	}

	//cout<<sqrt(get_mbb().diagonal_length())/HiMesh::sampled_points_num<<endl;;

	float dist = DBL_MAX;

	struct timeval start = get_cur_time();

	double smp_tm = 0;
	double usetree_tm = 0;
	double collect_triangle_tm = 0;
	double caldist_tm = 0;
	double ph_caldist_tm = 0.0;

	uint num_triangle = size_of_triangles();
	const float area_unit = area()/(num_triangle*sampling_rate);

	// associate each compressed facet with a list of original triangles, vice versa
	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){

		struct timeval start = get_cur_time();
		if(fit->rg==NULL || !fit->isSplittable()){
			continue;
		}
		fit->resetHausdorff();
		float fit_hdist = 0.0;
		fit->triangles = triangulate(fit);

		for(Triangle &cur_tri:fit->triangles) {
			float curhdist = 0.0;

			/* sampling the current triangle */
			unordered_set<Point> points;
			sample_points(cur_tri, points, area_unit);
			smp_tm += get_time_elapsed(start, true);

			// simply calculate against the BVH-tree built on the original mesh
			if(HiMesh::calculate_method == HCT_BVHTREE)
			{
				for(auto p:points){
					float dist = distance_tree(p);
					curhdist = max(curhdist, dist);
				}
				caldist_tm += get_time_elapsed(start, true);
			}

			// to collect the triangles associated
			// todo: this part of the latency is involved by an earlier
			// version of the implementation, will be removed in a future release
			// which will record all the removed facets instead of vertices
			vector<MyTriangle *> triangles;
			if(HiMesh::calculate_method != HCT_BVHTREE){
				for(const Point &p:fit->rg->removed_vertices){
					for(MyTriangle *t:VFmap[p]){
						if(!t->processed){
							t->processed = true;
							triangles.push_back(t);
						}
					}
				}
				for(const Point &p:fit->rg->removed_vertices){
					for(MyTriangle *t:VFmap[p]){
						t->processed = false;
					}
				}
				assert(triangles.size()>0);
				collect_triangle_tm += get_time_elapsed(start, true);
			}

			// use the triangles associated for calculating the Hasudorff distance
			if(HiMesh::calculate_method == HCT_ASSOCIATE)
			{
				// brute-forcely calculate
				for(auto p:points){
					float dist = DBL_MAX;
					for(MyTriangle *t:triangles){
						dist = min(dist, hispeed::PointTriangleDist((const float *)&p, (const float *)&t->tri));
					}
					curhdist = max(curhdist, dist);
				}
				caldist_tm += get_time_elapsed(start, true);
			}

			uint original_num = triangles.size();
			// further filter with the Triangle Cylinder
			if(HiMesh::calculate_method == HCT_ASSOCIATE_CYLINDER){
				for(vector<MyTriangle *>::iterator iter = triangles.begin();iter!=triangles.end();){
					bool keep = false;
					for(const Point &p:(*iter)->sampled_points){
						if(hispeed::PointInTriangleCylinder((const float *)&p, (const float *)&cur_tri)){
							keep = true;
							break;
						}
					}
					if(!keep){
						triangles.erase(iter);
					}else{
						iter++;
					}
				}
				collect_triangle_tm += get_time_elapsed(start, true);
			}

			for(MyTriangle *t:triangles){
				t->add_facet(fit);
			}

			if(HiMesh::calculate_method == HCT_ASSOCIATE_CYLINDER) {
				// brute-forcely calculate
				for(auto p:points){
					float dist = DBL_MAX;
					for(MyTriangle *t:triangles){
						dist = min(dist, hispeed::PointTriangleDist((const float *)&p, (const float *)&t->tri));
					}
					if(triangles.size() == 0){
						dist = 0;
					}
					curhdist = max(curhdist, dist);
				}
				caldist_tm += get_time_elapsed(start, true);
			}

			//log("#triangles: %ld-%d hdist: %.3f-%.3f-%.3f", original_num, triangles.size(), curhdist[0], curhdist[1], curhdist[2]);

			// collect results
			fit_hdist = max(fit_hdist, curhdist);
			triangles.clear();
		}
		// update the hausdorff distance
		fit->updateHausdorff(fit_hdist);
	}
	start = get_cur_time();
	/*
	 *
	 * for computing the proxy hausdorf distance
	 *
	 * */
	if(HiMesh::calculate_method == HCT_BVHTREE){
		list<Triangle> triangles = get_triangles();
		TriangleTree *tree = new TriangleTree(triangles.begin(), triangles.end());
		tree->build();
		tree->accelerate_distance_queries();
		float phdist = 0.0;
		struct timeval ss = get_cur_time();
		for(MyTriangle *t:original_facets){
			// for each sampled point, find the closest facet to it
			// and it will be proxy facet of that point
			for(const Point &p:t->sampled_points){
				float dist = tree->squared_distance(p);
				phdist = max(phdist, dist);
			}
		}
		logt("BVH %f", ss, sqrt(phdist));
		phdist = sqrt(phdist);
	}

	if(HiMesh::calculate_method == HCT_ASSOCIATE || HiMesh::calculate_method == HCT_ASSOCIATE_CYLINDER){
		for(MyTriangle *t:original_facets){
			tt++;
			//log("%d", t->facets.size());
			if(t->facets.size()>0){
				Face_iterator cur_fit;
				Point cur_p;
				// calculate the Hausdorf distance
				float hdist = 0.0;
				// for each sampled point, find the closest facet to it
				// and it will be proxy facet of that point
				for(const Point &p:t->sampled_points){
					float dist = DBL_MAX;
					Face_iterator fit;
					for(Face_iterator &ft:t->facets){
						for(Triangle &tt:ft->triangles){
							float d = hispeed::PointTriangleDist((const float *)&p, (const float *)&tt);
							if(d<dist){
								dist = d;
								fit = ft;
							}
						}
					}
					fit->updateProxyHausdorff(dist);
					// get the maximum
					if(hdist < dist){
						hdist = dist;
						cur_p = p;
						cur_fit = fit;
					}
				}
//
//				//if(t->facets.size()>1)
//				if(false && (tt == 32149 || t->id == 2)){
//					vector<Triangle> compressed;
//					for(HiMesh::Face_iterator fit:t->facets){
//						compressed.insert(compressed.end(), fit->triangles.begin(), fit->triangles.end());
//					}
//					hispeed::write_triangles(compressed, "/gisdata/a.compressed.off");
//
//					vector<Triangle> origin;
//					origin.push_back(t->tri);
//					hispeed::write_triangles(origin, "/gisdata/a.origin.off");
//
//					hispeed::write_triangles(cur_fit->triangles, "/gisdata/a.selected.off");
//
//					vector<Point> points;
//					points.push_back(cur_p);
//					hispeed::write_points(points, "/gisdata/a.points.off");
//
//					vector<Point> sampledpoints;
//					sampledpoints.assign(t->sampled_points.begin(), t->sampled_points.end());
//					hispeed::write_points(sampledpoints, "/gisdata/a.sampled.off");
//
//					write_to_off("/gisdata/a.current.off");
//					log("%d %f %d", tt, hdist, t->id);
//					exit(0);
//				}
			}
			t->reset();
		}
	}
	ph_caldist_tm += get_time_elapsed(start, true);

	pair<float, float> current_hausdorf = pair<float, float>(0.0, 0.0);
	// some collected statistics for presentation
	float min_hdist = DBL_MAX;
	float max_hdist = 0.0;
	float avg_hdist = 0.0;
	float min_proxy_hdist = DBL_MAX;
	float max_proxy_hdist = 0.0;
	float avg_proxy_hdist = 0.0;
	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
		const float fit_hdist = fit->getHausdorff()+sqrt(area_unit)/sqrt(2);
		const float fit_proxy_hdist = fit->getProxyHausdorff()+sqrt(area_unit)/sqrt(2);

		fit->setHausdorff(fit_hdist);
		fit->setProxyHausdorff(fit_proxy_hdist);

		min_hdist = min(min_hdist, fit_hdist);
		max_hdist = max(max_hdist, fit_hdist);
		avg_hdist += fit_hdist;
		min_proxy_hdist = min(min_proxy_hdist, fit_proxy_hdist);
		max_proxy_hdist = max(max_proxy_hdist, fit_proxy_hdist);
		avg_proxy_hdist += fit_proxy_hdist;
		current_hausdorf.second = max(current_hausdorf.second, fit->getHausdorff());
		current_hausdorf.first = max(current_hausdorf.first, fit->getProxyHausdorff());
//		if(this->size_of_facets()<100){
//			log("%f\t%f", fit->getProxyHausdorff(), fit->getHausdorff());
//		}
	}
	avg_hdist /= (int)size_of_facets();
	avg_proxy_hdist /= (int)size_of_facets();
	globalHausdorfDistance.push_back(current_hausdorf);

	if(global_ctx.verbose>=2)
	{
		log("step: %2d smp: %.3f tri: %.3f h_cal: %.3f ph_cal:%.3f avg_hdist#vertices: %ld #facets: %ld",
				i_curDecimationId,
				smp_tm, collect_triangle_tm, caldist_tm,ph_caldist_tm,
				size_of_vertices(), size_of_triangles());
		log("hausdorff(min-avg-max): %.5f-%.5f-%.5f	proxy_hausdorff(min-avg-max): %.5f-%.5f-%.5f	sampling_rate: %.10f",
				min_hdist, avg_hdist, max_hdist,
				min_proxy_hdist, avg_proxy_hdist, max_proxy_hdist,
				sqrt(area_unit)/sqrt(2));
		log("encode %d:\t[%.2f %.2f]\t%ld", i_curDecimationId, current_hausdorf.first, current_hausdorf.second, size_of_vertices());
	}
}


}
