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
#include "../PPMC/configuration.h"
#include "../PPMC/mymesh.h"
#include "geometry.h"

float point_to_face_distance(Point p, MyMesh::Face_iterator fit){
	float triangle[9];
	float point[3];
	point[0] = p.x();
	point[1] = p.y();
	point[2] = p.z();
	MyMesh::Halfedge_const_handle hd = fit->halfedge();
	MyMesh::Halfedge_const_handle h = hd->next();
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

pair<float, float> MyMesh::getHausdorfDistance(){
	assert(i_nbDecimations>=i_curDecimationId);
	return i_nbDecimations>i_curDecimationId?globalHausdorfDistance[i_nbDecimations - i_curDecimationId-1]:std::pair<float, float>(0, 0);
}

pair<float, float> MyMesh::getNextHausdorfDistance(){
	assert(i_nbDecimations>i_curDecimationId);
	return i_nbDecimations>(i_curDecimationId+1)?globalHausdorfDistance[i_nbDecimations - i_curDecimationId - 2]:std::pair<float, float>(0, 0);
}

/**
  * Start the next compression operation.
  */
void MyMesh::startNextCompresssionOp()
{
	for(MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
		vit->resetState();

	for(MyMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
		hit->resetState();

	for(MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
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
void MyMesh::decimationStep()
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
MyMesh::Halfedge_handle MyMesh::vertexCut(Halfedge_handle startH)
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

void MyMesh::computeHausdorfDistance(){

	pair<float, float> current_hausdorf = pair<float, float>(0.0, 0.0);
	float dist = DBL_MAX;

	for(MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
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
		current_hausdorf.second = max(current_hausdorf.second, maxdist);
		dist = min(dist, maxdist);
	}
	globalHausdorfDistance.push_back(current_hausdorf);
	if(global_ctx.verbose>=2){
		log("encode %d:\t[%.2f %.2f]\t%ld", i_curDecimationId, dist, current_hausdorf.second, size_of_vertices());
	}
}

/**
  * One step of the removed vertex coding conquest.
  */
void MyMesh::RemovedVertexCodingStep()
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
void MyMesh::InsertedEdgeCodingStep()
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
void MyMesh::encodeInsertedEdges(unsigned i_operationId)
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
void MyMesh::encodeRemovedVertices(unsigned i_operationId)
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
