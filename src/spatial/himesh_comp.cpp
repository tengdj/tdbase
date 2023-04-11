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


void HiMesh::encode(int lod){
	assert(is_compression_mode());
	b_jobCompleted = false;
	while(!b_jobCompleted) {
		startNextCompresssionOp();
	}
}

/**
  * Start the next compression operation.
  */
void HiMesh::startNextCompresssionOp()
{
	// 1. reset the stats
	for(HiMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
		vit->resetState();

	for(HiMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
		hit->resetState();

	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
		fit->resetState();
	// Reset the number of removed vertices.
	i_nbRemovedVertices = 0;


	// 2. do one round of decimation
	//choose a halfedge that can be processed:
	{
		// teng: we always start from the middle, DO NOT use the rand function
		// size_t i_heInitId = (float)rand() / RAND_MAX * size_of_halfedges();
		size_t i_heInitId = size_of_halfedges()/2;
		Halfedge_iterator hitInit = halfedges_begin();
//		for (unsigned i = 0; i < i_heInitId; ++i)
//		  ++hitInit;
		hitInit->setInQueue();
		gateQueue.push((Halfedge_handle)hitInit);
	}
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

	// 3. do the encoding job
    if (i_nbRemovedVertices == 0){
		b_jobCompleted = true;
		i_nbDecimations = i_curDecimationId--;

		// Write the compressed data to the buffer.
	    writeBaseMesh();
	    int i_deci = i_curDecimationId;
	    assert(i_deci>0);
	    while (i_deci>=0)
	    {
			encodeRemovedVertices(i_deci);
			encodeInsertedEdges(i_deci);
			i_deci--;
	    }
	} else {
		// 3dpro: compute the hausdorf distance for all the existing faces
		computeHausdorfDistance();

		RemovedVertexCodingStep();
		InsertedEdgeCodingStep();
		// finish this round of decimation and start the next
	    i_curDecimationId++; // Increment the current decimation operation id.
	}
}

HiMesh::replacing_group *merge(vector<HiMesh::replacing_group *> &reps){

	HiMesh::replacing_group *ret = reps[0];
	for(size_t i=1;i<reps.size();i++){
		for(auto a:reps[i]->added_faces){
			ret->added_faces.emplace(a);
		}
		for(auto a:reps[i]->removed_vertices){
			ret->removed_vertices.emplace(a);
		}
		delete reps[i];
	}
	reps.clear();
	return ret;
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
	vector<Face_handle> new_faces;
	vector<replacing_group *> rep_groups;
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

		  Face_handle new_face = hCorner->face();
		  Face_handle old_face = hCorner->opposite()->face();

		  new_faces.push_back(new_face);

		  // the old face will be removed in the vertex cut
		  if(map_group.find(old_face)!=map_group.end()){
			  rep_groups.push_back(map_group[old_face]);
			  map_group[old_face]->added_faces.erase(old_face);
			  map_group.erase(old_face);
		  }
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
	Face_handle added_face = hNewFace->facet();
	new_faces.push_back(added_face);

	log("%ld %ld", rep_groups.size(), new_faces.size());
	replacing_group *rg = NULL;
	if(rep_groups.size()==0){
		rg = new replacing_group();
	}else{
		// the faces are
		rg = merge(rep_groups);
	}
	assert(rg);
	for(Face_handle nf:new_faces){
		rg->added_faces.emplace(nf);
		map_group[nf] = rg;
	}
	rg->removed_vertices.emplace(vPos);

	//now mark the new face as having a removed vertex
	added_face->setSplittable();
	// keep the removed vertex position.
	added_face->setRemovedVertexPos(vPos);
	added_face->addImpactPoints(impactpoints);
	added_face->addImpactPoint(vPos);
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
        	// record the removed points during compressing.
            removedPoints.push_back(rmved);
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
	// Resize the vector to add the current conquest symbols.
	connectEdgeSym.push_back(std::deque<unsigned>());
	// Add the first halfedge to the queue.
	pushHehInit();
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

// Write the base mesh.
void HiMesh::writeBaseMesh()
{
    // Write the bounding box min coordinate.
    for (unsigned i = 0; i < 3; ++i)
        writeFloat((float)(mbb.low[i]));
    for (unsigned i = 0; i < 3; ++i)
        writeFloat((float)(mbb.high[i]));
    // Write the quantization step.

    unsigned i_nbVerticesBaseMesh = size_of_vertices();
    unsigned i_nbFacesBaseMesh = size_of_facets();

    // Write the number of level of decimations.
    writeInt16(i_nbDecimations);

    // Write the number of vertices and faces on 16 bits.
    writeInt(i_nbVerticesBaseMesh);
    writeInt(i_nbFacesBaseMesh);

    // Write the vertices of the edge that is the departure of the coding conquests.
    size_t id = 0;
    for (unsigned j = 0; j < 2; ++j) {
    	writePoint(vh_departureConquest[j]->point());
        vh_departureConquest[j]->setId(id++);
    }

    // Write the other vertices.
    for (HiMesh::Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    {
        if (vit == vh_departureConquest[0] || vit == vh_departureConquest[1])
            continue;
        writePoint(vit->point());
        // Set an id to the vertex.
        vit->setId(id++);
    }

    // Write the base mesh face vertex indices.
    for (HiMesh::Facet_iterator fit = facets_begin();
         fit != facets_end(); ++fit)
    {
        unsigned i_faceDegree = fit->facet_degree();
        writeInt(i_faceDegree);
        Halfedge_around_facet_const_circulator hit(fit->facet_begin()), end(hit);
        do
        {
            // Write the current vertex id.
        	writeInt(hit->vertex()->getId());
        }
        while(++hit != end);
    }

    // Write the hausdorf distance for the base faces
    for (HiMesh::Facet_iterator fit = facets_begin();
         fit != facets_end(); ++fit)
    {
        writeFloat(fit->getHausdorfDistance().first);
        writeFloat(fit->getHausdorfDistance().second);
    }

    // 3dpro
    // Write the maximum volume change for each round of decimation
    assert(globalHausdorfDistance.size()==i_nbDecimations);
    for(unsigned i=0;i<i_nbDecimations;i++){
    	writeFloat(globalHausdorfDistance[i].first);
    	writeFloat(globalHausdorfDistance[i].second);
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
    for (unsigned i = 0; i < i_lenConn; ++i) {
        // Encode the connectivity.
        unsigned sym = connSym[i];
        writeChar(sym);
        // Encode the geometry if necessary.
        if (sym) {
            writePoint(geomSym[k]);
            k++;
        }
        writeuInt16(hausSym[i]);
    }
}

}
