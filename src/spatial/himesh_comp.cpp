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

namespace tdbase{


void HiMesh::encode(){
	assert(is_compression_mode());
	b_jobCompleted = false;
	while(!b_jobCompleted) {
		startNextCompresssionOp();
	}
}

/**
  * Push the first halfedge for the coding and decoding conquest in the gate queue.
  */
void HiMesh::pushHehInit() {
    // Find the first halfedge.
    Halfedge_handle hehBegin;
    Halfedge_around_vertex_circulator hit(vh_departureConquest[0]->vertex_begin());
    while (true) {
        hehBegin = hit->opposite();
        if (hehBegin->vertex() == vh_departureConquest[1])
            break;
        ++hit;
    }
    // Push it to the queue.
    gateQueue.push(hehBegin);
}


/**
  * Test if a vertex is removable.
  */
bool HiMesh::isRemovable(Vertex_handle v) const {
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
	  } while(++hit != end);
	  //
	  bool removable = !willViolateManifold(heh_oneRing);
	  //removable &= isConvex(vh_oneRing);
	  return removable;
	}
	return false;
}

/**
  * Start the next compression operation.
  */
void HiMesh::startNextCompresssionOp() {
	// 1. reset the stats
	for(HiMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
		vit->resetState();

	for(HiMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
		hit->resetState();

	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
		fit->resetState();

	i_nbRemovedVertices = 0; // Reset the number of removed vertices.

	// 2. do one round of decimation
	//choose a halfedge that can be processed:
	if(i_curDecimationId<10) {
		// teng: we always start from the middle, DO NOT use the rand function
		// size_t i_heInitId = (float)rand() / RAND_MAX * size_of_halfedges();
		size_t i_heInitId = size_of_halfedges()/2;
		Halfedge_iterator hitInit = halfedges_begin();
		for (unsigned i = 0; i < i_heInitId; ++i)
		  ++hitInit;
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
		encodeBaseMesh();
	    int i_deci = i_curDecimationId;
	    assert(i_deci>0);
	    while (i_deci>=0) {
			encodeHausdorff(i_deci);
			encodeRemovedVertices(i_deci);
			encodeInsertedEdges(i_deci);
			i_deci--;
	    }
	} else {
		// tdbase: compute and encode the Hausdorff distance for all the facets in this LOD
		computeHausdorfDistance();
		HausdorffCodingStep();
		// finish this round of decimation and start the next
		RemovedVertexCodingStep();
		InsertedEdgeCodingStep();
	    i_curDecimationId++; // Increment the current decimation operation id.
	}
}

/**
  * Perform the re-edging and the vertex removal.
  * Store the position of the removed vertex.
  * TODO: critical
  */
HiMesh::Halfedge_handle HiMesh::vertexCut(Halfedge_handle startH) {
	Vertex_handle v = startH->vertex();

	//make sure that the center vertex can be removed
	assert(!v->isConquered());
	assert(v->vertex_degree()>2);

	Halfedge_handle h = startH->opposite(), end(h);
	int removed = 0;
	do
	{
		assert(!h->is_border());
		Face_handle f = h->facet();
		assert(!f->isConquered()); //we cannot cut again an already cut face, or a NULL patch

		//if the face is not a triangle, cut the corner to make it a triangle
		if(f->facet_degree()>3)	{
			//loop around the face to find the appropriate other halfedge
			Halfedge_handle hSplit(h->next());
			for(; hSplit->next()->next() != h; hSplit = hSplit->next())
				;
			Halfedge_handle hCorner = split_facet(h, hSplit);
			//mark the new halfedges as added
			hCorner->setAdded();
			hCorner->opposite()->setAdded();
			// the corner one inherit the original facet
			// while the fRest is a newly generated facet
			Face_handle fCorner = hCorner->face();
			Face_handle fRest = hCorner->opposite()->face();
			//log("split %ld + %ld %ld", fCorner->facet_degree(), fRest->facet_degree(), f->facet_degree());
		}
		//mark the vertex as conquered
		h->vertex()->setConquered();
		removed++;
	} while((h=h->opposite()->next()) != end);

	//copy the position of the center vertex:
	Point vPos = startH->vertex()->point();

	int bf = size_of_facets();
	//remove the center vertex
	Halfedge_handle hNewFace = erase_center_vertex(startH);
	Face_handle added_face = hNewFace->facet();

	//log("test: %d = %d - %ld merged %ld replacing groups", removed, bf, size_of_facets(), rep_groups.size());

	//now mark the new face as having a removed vertex
	added_face->setSplittable();
	// keep the removed vertex position.
	added_face->setRemovedVertexPos(vPos);

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
void HiMesh::RemovedVertexCodingStep() {
    // Resize the vectors to add the current conquest symbols.
    geometrySym.push_back(std::deque<Point>());
    connectFaceSym.push_back(std::deque<unsigned>());

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
void HiMesh::InsertedEdgeCodingStep() {
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

void HiMesh::HausdorffCodingStep(){

	// we store two times for experimental evaluation
	hausdorfSym.push_back(std::deque<unsigned char>());
	proxyhausdorfSym.push_back(std::deque<unsigned char>());
	hausdorfSym_float.push_back(std::deque<float>());
	proxyhausdorfSym_float.push_back(std::deque<float>());
	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){

		// Determine the hausdorf symbol
		// the hausdorf distances for this round of decimation.
		const pair<float, float> current_hausdorf = globalHausdorfDistance[globalHausdorfDistance.size()-1];
		// 3dpro, besides the symbol, we also encode the hausdorf distance into it.
		assert(fit->getProxyHausdorff()<=current_hausdorf.first);
		assert(fit->getHausdorff()<=current_hausdorf.second);

		unsigned char proxyhausdorf = current_hausdorf.first==0?0:ceil(fit->getProxyHausdorff()/current_hausdorf.first*127);
		unsigned char hausdorf = current_hausdorf.second==0?0:ceil(fit->getHausdorff()/current_hausdorf.second*127);
		hausdorfSym[i_curDecimationId].push_back(hausdorf);
		proxyhausdorfSym[i_curDecimationId].push_back(proxyhausdorf);
		hausdorfSym_float[i_curDecimationId].push_back(fit->getHausdorff());
		proxyhausdorfSym_float[i_curDecimationId].push_back(fit->getProxyHausdorff());
		if(global_ctx.verbose>=3)
		{
			log("encode facet: %.2f %.2f %d %d", fit->getProxyHausdorff(), fit->getHausdorff(), proxyhausdorf, hausdorf);
		}
	}
}

// Write the base mesh.
void HiMesh::encodeBaseMesh() {
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

    // 3dpro
    // Write the mesh-level Hausdorf
    assert(globalHausdorfDistance.size()==i_nbDecimations);
    for(int i=i_nbDecimations-1;i>=0;i--){
    	writeFloat(globalHausdorfDistance[i].first);
    	writeFloat(globalHausdorfDistance[i].second);
    	//log("encode hausdorff: %d %f %f", i, globalHausdorfDistance[i].first, globalHausdorfDistance[i].second);
    }
}

/**
  * Encode an inserted edge list.
  */
void HiMesh::encodeInsertedEdges(unsigned i_operationId) {

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
void HiMesh::encodeRemovedVertices(unsigned i_operationId) {
    std::deque<unsigned> &connSym = connectFaceSym[i_operationId];
    std::deque<Point> &geomSym = geometrySym[i_operationId];

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
    }
}

void HiMesh::encodeHausdorff(unsigned i_operationId){
	std::deque<unsigned char> &hausSym = hausdorfSym[i_operationId];
	std::deque<unsigned char> &proxyhausSym = proxyhausdorfSym[i_operationId];
	std::deque<float> &hausSym_float = hausdorfSym_float[i_operationId];
	std::deque<float> &proxyhausSym_float = proxyhausdorfSym_float[i_operationId];
	for(int i=0;i<hausSym.size();i++){
		writeChar(hausSym[i]);
		writeChar(proxyhausSym[i]);
		writeFloat(hausSym_float[i]);
		writeFloat(proxyhausSym_float[i]);
	}
}

}
