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
#include <CGAL/squared_distance_3.h>

#include "math.h"
#include "himesh.h"
#include "geometry.h"

namespace hispeed{


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

#define PI 3.1415926
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
			double angle = acos(cosvalue)*180/PI;
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


float point_to_face_distance(Point p, HiMesh::Face_iterator fit){
	float triangle[9];
	float point[3];
	point[0] = p.x();
	point[1] = p.y();
	point[2] = p.z();
	HiMesh::Halfedge_const_handle hd = fit->halfedge();
	HiMesh::Halfedge_const_handle h = hd->next();
	float mindist = DBL_MAX;
	while(h->next()!=hd){
		Point p1 = hd->vertex()->point();
		Point p2 = h->vertex()->point();
		Point p3 = h->next()->vertex()->point();
		h = h->next();
		triangle[0] = p1.x();
		triangle[1] = p1.y();
		triangle[2] = p1.z();
		triangle[3] = p2.x();
		triangle[4] = p2.y();
		triangle[5] = p2.z();
		triangle[6] = p3.x();
		triangle[7] = p3.y();
		triangle[8] = p3.z();
		float dist = PointTriangleDist((const float*)point, (const float*)triangle);
		mindist = min(mindist, dist);
	}
	return mindist;
}

pair<float, float> HiMesh::getHausdorfDistance(){
	assert(i_nbDecimations>=i_curDecimationId);
	return i_nbDecimations>i_curDecimationId?globalHausdorfDistance[i_nbDecimations - i_curDecimationId-1]:std::pair<float, float>(0, 0);
}

pair<float, float> HiMesh::getNextHausdorfDistance(){
	assert(i_nbDecimations>i_curDecimationId);
	return i_nbDecimations>(i_curDecimationId+1)?globalHausdorfDistance[i_nbDecimations - i_curDecimationId - 2]:std::pair<float, float>(0, 0);
}

/**
  * Start the next compression operation.
  */
void HiMesh::startNextCompresssionOp()
{
	for(HiMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
		vit->resetState();

	for(HiMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
		hit->resetState();

	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
		fit->resetState();

    decimationStep();
    if (i_nbRemovedVertices == 0){
		b_jobCompleted = true;
		i_curDecimationId--;
		writeCompressedData();
	} else {
		// 3dpro: compute the hausdorf distance for all the existing faces
		computeHausdorfDistance();

		RemovedVertexCodingStep();
		InsertedEdgeCodingStep();
		// finish this round of decimation and start the next
	    i_curDecimationId++; // Increment the current decimation operation id.
	}
}

// One decimation step.
void HiMesh::decimationStep()
{
	// Select the first gate to begin the decimation.
	// teng: we always start from the middle, DO NOT use the pushHehInit() function
	// size_t i_heInitId = (float)rand() / RAND_MAX * size_of_halfedges();
	size_t i_heInitId = size_of_halfedges()/2;
	Halfedge_iterator hitInit = halfedges_begin();
	for (unsigned i = 0; i < i_heInitId; ++i)
	  ++hitInit;
	hitInit->setInQueue();
	gateQueue.push((Halfedge_handle)hitInit);

	// Reset the number of removed vertices.
	i_nbRemovedVertices = 0;

    //choose a halfedge that can be processed:
    while(!gateQueue.empty()) {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        //pick a the first halfedge from the queue. f is the adjacent face.
        assert(!h->is_border());
        Face_handle f = h->facet();

        //if the face is already processed, pick the next halfedge:
        if(f->isConquered())
        {
            h->removeFromQueue();
            continue;
        }
        //the face is not processed. Count the number of non conquered vertices that can be split
        bool hasRemovable = false;
        Halfedge_handle unconqueredVertexHE;

        for(Halfedge_handle hh = h->next(); hh!=h; hh=hh->next())
        {
            if(isRemovable(hh->vertex()))
            {
            	hasRemovable = true;
                unconqueredVertexHE = hh;
                break;
            }
        }

        //if all face vertices are conquered, then the current face is a null patch:
        if(!hasRemovable)
        {
            f->setUnsplittable();
            //and add the outer halfedges to the queue. Also mark the vertices of the face conquered
            Halfedge_handle hh = h;
            do
            {
                hh->vertex()->setConquered();
                Halfedge_handle hOpp = hh->opposite();
                assert(!hOpp->is_border());
                if(!hOpp->facet()->isConquered())
                {
                    gateQueue.push(hOpp);
                    hOpp->setInQueue();
                }
            }
            while((hh = hh->next()) != h);
            h->removeFromQueue();
        } else {
            //in that case, cornerCut that vertex.
            h->removeFromQueue();
            vertexCut(unconqueredVertexHE);
        }
    }
}


/**
  * Perform the re-edging and the vertex removal.
  * Store the position of the removed vertex.
  */
HiMesh::Halfedge_handle HiMesh::vertexCut(Halfedge_handle startH)
{
	Vertex_handle v = startH->vertex();

	//make sure that the center vertex can be removed
	assert(!v->isConquered());
	assert(v->vertex_degree()>2);

	vector<Point> impactpoints;
	//int i = 0;
	Halfedge_handle h = startH->opposite(), end(h);
	do
	{
		assert(!h->is_border());
		Face_handle f = h->facet();
		assert(!f->isConquered()); //we cannot cut again an already cut face, or a NULL patch

		vector<Point> ips = f->getImpactPoints();
		impactpoints.insert(impactpoints.end(), ips.begin(), ips.end());

		//if the face is not a triangle, cut the corner
		int deg_bef = f->facet_degree();
		if(f->facet_degree()>3)
		{
		  //loop around the face to find the appropriate other halfedge
		  Halfedge_handle hSplit(h->next());
		  for(; hSplit->next()->next() != h; hSplit = hSplit->next())
				;
		  Halfedge_handle hCorner = split_facet(h, hSplit);
		  //mark the new halfedges as added
		  hCorner->setAdded();
		  hCorner->opposite()->setAdded();

		  //hCorner->opposite()->facet()->resetImpactPoints();
		  //log("addr: %ld %ld %ld", f, hCorner->facet(), hCorner->opposite()->facet());
		  //log("split: %ld %ld", hCorner->facet()->getImpactPoints().size(), hCorner->opposite()->facet()->getImpactPoints().size());
		  //hCorner->facet()->addImpactPoints(ips);

		  //log("%d %d",deg_bef, hCorner->opposite()->facet_degree());
		}
		//mark the vertex as conquered
		h->vertex()->setConquered();
	} while((h=h->opposite()->next()) != end);

	//copy the position of the center vertex:
	Point vPos = startH->vertex()->point();

	//remove the center vertex
	Halfedge_handle hNewFace = erase_center_vertex(startH);

	//now mark the new face as having a removed vertex
	hNewFace->facet()->setSplittable();
	// keep the removed vertex position.
	hNewFace->facet()->setRemovedVertexPos(vPos);
	hNewFace->facet()->addImpactPoints(impactpoints);
	hNewFace->facet()->addImpactPoint(vPos);
	//scan the outside halfedges of the new face and add them to
	//the queue if the state of its face is unknown. Also mark it as in_queue
	h = hNewFace;
	do
	{
		Halfedge_handle hOpp = h->opposite();
		assert(!hOpp->is_border());
		if(!hOpp->facet()->isConquered())
		{
			gateQueue.push(hOpp);
			hOpp->setInQueue();
		}
	}
	while((h = h->next()) != hNewFace);

	// Increment the number of removed vertices.
	i_nbRemovedVertices++;

	return hNewFace;
}

void HiMesh::computeHausdorfDistance(){

	pair<float, float> current_hausdorf = pair<float, float>(0.0, 0.0);
	float dist = DBL_MAX;

	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
		if(fit->isSplittable()){
			vector<Point> ps = fit->getImpactPoints();
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
		}
		current_hausdorf.second = max(current_hausdorf.second, fit->getHausdorfDistance().second);
		dist = min(dist, fit->getHausdorfDistance().second);
	}
	globalHausdorfDistance.push_back(current_hausdorf);
	if(global_ctx.verbose>=2){
		log("encode %d:\t[%.2f %.2f]\t%ld", i_curDecimationId, dist, current_hausdorf.second, size_of_vertices());
	}
}

/**
  * One step of the removed vertex coding conquest.
  */
void HiMesh::RemovedVertexCodingStep()
{
    // Resize the vectors to add the current conquest symbols.
    geometrySym.push_back(std::deque<Point>());
    connectFaceSym.push_back(std::deque<unsigned>());
    hausdorfSym.push_back(std::deque<unsigned>());

    // Add the first halfedge to the queue.
    pushHehInit();
    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        Face_handle f = h->facet();

        // If the face is already processed, pick the next halfedge:
        if (f->isProcessed())
            continue;

        // Determine face symbol.
        unsigned sym = f->isSplittable();

        // Push the symbols.
        connectFaceSym[i_curDecimationId].push_back(sym);

        // Determine the geometry symbol.
        if (sym){
            Point rmved = f->getRemovedVertexPos();
            geometrySym[i_curDecimationId].push_back(rmved);
        }

        // Determine the hausdorf symbol
        // the hausdorf distances for this round of decimation.
		pair<float, float> current_hausdorf = this->globalHausdorfDistance[this->globalHausdorfDistance.size()-1];
		// 3dpro, besides the symbol, we also encode the hausdorf distance into it.
		unsigned con = current_hausdorf.first==0?0:(f->getHausdorfDistance().first/current_hausdorf.first*100.0);
		unsigned pro = current_hausdorf.second==0?0:(f->getHausdorfDistance().second/current_hausdorf.second*100.0);
        hausdorfSym[i_curDecimationId].push_back((con<<8)|pro);
        if(global_ctx.verbose>=3){
        	log("encode face: %d %.2f %d %.2f %d",sym, f->getHausdorfDistance().first, con, f->getHausdorfDistance().second, pro);
        }

        // Mark the face as processed.
        f->setProcessedFlag();

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h;
        do
        {
            Halfedge_handle hOpp = hIt->opposite();
            assert(!hOpp->is_border());
            if (!hOpp->facet()->isProcessed())
                gateQueue.push(hOpp);
            hIt = hIt->next();
        }
        while (hIt != h);
    }
}

/**
  * One step of the inserted edge coding conquest.
  */
void HiMesh::InsertedEdgeCodingStep()
{

	// Add the first halfedge to the queue.
	pushHehInit();
	// Resize the vector to add the current conquest symbols.
	connectEdgeSym.push_back(std::deque<unsigned>());
    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        // Test if the edge has already been conquered.
        if (h->isProcessed())
            continue;

        // Mark the halfedge as processed.
        h->setProcessed();
        h->opposite()->setProcessed();

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h->next();
        while (hIt->opposite() != h)
        {
            if (!hIt->isProcessed())
                gateQueue.push(hIt);
            hIt = hIt->opposite()->next();
        }

        // Don't write a symbol if the two faces of an edgde are unsplitable.
        // this can help to save some space, since it is guaranteed that the edge is not inserted
        bool b_toCode = h->facet()->isUnsplittable()
                        && h->opposite()->facet()->isUnsplittable()
                        ? false : true;

        // Determine the edge symbol.
        unsigned sym;
        if (h->isOriginal())
            sym = 0;
        else
            sym = 1;

        // Store the symbol if needed.
        if (b_toCode)
            connectEdgeSym[i_curDecimationId].push_back(sym);
    }
}


/**
  * Encode an inserted edge list.
  */
void HiMesh::encodeInsertedEdges(unsigned i_operationId)
{

    std::deque<unsigned> &symbols = connectEdgeSym[i_operationId];
    assert(symbols.size() > 0);

    unsigned i_len = symbols.size();
    for (unsigned i = 0; i < i_len; ++i)
    {
        writeChar(symbols[i]);
    }
}


/**
  * Encode the geometry and the connectivity of a removed vertex list.
  */
void HiMesh::encodeRemovedVertices(unsigned i_operationId)
{
    std::deque<unsigned> &connSym = connectFaceSym[i_operationId];
    std::deque<Point> &geomSym = geometrySym[i_operationId];
    std::deque<unsigned> &hausSym = hausdorfSym[i_operationId];

    unsigned i_lenGeom = geomSym.size();
    unsigned i_lenConn = connSym.size();
    assert(i_lenGeom > 0);
    assert(i_lenConn > 0);

    unsigned k = 0;
    for (unsigned i = 0; i < i_lenConn; ++i)
    {
        // Encode the connectivity.
        unsigned sym = connSym[i];
        writeChar(sym);
        // Encode the geometry if necessary.
        if (sym)
        {
            Point p = geomSym[k];
            for (unsigned j = 0; j < 3; ++j)
            {
            	writeFloat(p[j]);
            }
            k++;
        }
        writeuInt16(hausSym[i]);
    }
}

}
