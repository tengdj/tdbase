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

#include "himesh.h"

namespace hispeed{

/**
  * Start the next decompression operation.
  */
void HiMesh::startNextDecompresssionOp()
{
	// reset the states
    for (HiMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
        hit->resetState();

    for (HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
        fit->resetState();

    if(global_ctx.verbose>=2 && ((float)i_curDecimationId / i_nbDecimations * 100 < i_decompPercentage||i_curDecimationId == i_nbDecimations)){
    	log("decode %d:\t%.2f\%\t[%.2f, %.2f]", i_curDecimationId, (float)i_curDecimationId / i_nbDecimations * 100, getHausdorfDistance().first, getHausdorfDistance().second);
    }
    if ((float)i_curDecimationId / i_nbDecimations * 100 >= i_decompPercentage){
    	if(i_curDecimationId == i_nbDecimations){
    		// reset all the hausdorf distance to 0 for the highest LOD
    		// as we do not have another round of decoding to set them
            for (HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
            	fit->setProgressive(0.0);
            	fit->setConservative(0.0);
            }
    	}
        b_jobCompleted = true;
    } else {

        undecimationStep();
		InsertedEdgeDecodingStep();
		insertRemovedVertices();
		removeInsertedEdges();

		i_curDecimationId++; // increment the current decimation operation id.
    }
}

/**
  * One undecimation step.
  */
void HiMesh::undecimationStep()
{
	// Add the first halfedge to the queue.
	pushHehInit();
	while (!gateQueue.empty())
	{
		Halfedge_handle h = gateQueue.front();
		gateQueue.pop();

		Face_handle f = h->facet();

		// If the face is already processed, pick the next halfedge:
		if (f->isConquered())
			continue;

		// Decode the face symbol.
		unsigned sym = readChar();

		// Add the other halfedges to the queue
		Halfedge_handle hIt = h;
		do
		{
			Halfedge_handle hOpp = hIt->opposite();
			assert(!hOpp->is_border());
			if (!hOpp->facet()->isConquered())
				gateQueue.push(hOpp);
			hIt = hIt->next();
		} while (hIt != h);
		// Decode the geometry symbol.
		if (sym == 1){
			float coord[3];
			for (unsigned i = 0; i < 3; ++i)
			{
				// Store the value.
				coord[i] = readFloat();
			}

			Point rmved(coord[0], coord[1], coord[2]);
			f->setSplittable();
			f->setRemovedVertexPos(rmved);
		} else {
			f->setUnsplittable();
		}

		// decode the hausdorf distance symbols
		uint hsym = readuInt16();
		uint prosym = hsym&((1<<8)-1);
		hsym >>= 8;
		uint consym = hsym&((1<<8)-1);
		pair<float, float> hdist = getNextHausdorfDistance();
		float progressive = prosym/100.0 * hdist.second;
		float conservative = consym/100.0 * hdist.first;
		f->setProgressive(progressive);
		f->setConservative(conservative);
		if(global_ctx.verbose>=3){
			log("decode face: %d %.2f %d %.2f %d", sym, conservative, consym, progressive, prosym);
		}
	}
}

/**
  * One step of the inserted edge coding conquest.
  */
void HiMesh::InsertedEdgeDecodingStep()
{
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

        // Test if there is a symbol for this edge.
        // There is no symbol if the two faces of an edge are unsplitable.
        if (h->facet()->isSplittable() || h->opposite()->facet()->isSplittable())
        {
            // Decode the edge symbol.
            unsigned sym = readChar();
            // Determine if the edge is original or not.
            // Mark the edge to be removed.
            if (sym != 0)
                h->setAdded();
        }

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h->next();
        while (hIt->opposite() != h)
        {
            if (!hIt->isProcessed() && !hIt->isNew())
                gateQueue.push(hIt);
            hIt = hIt->opposite()->next();
        }
        assert(!hIt->isNew());

    }
}


/**
  * Insert center vertices.
  */
void HiMesh::insertRemovedVertices()
{

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

        assert(!h->isNew());

        if (f->isSplittable())
        {
            // Insert the vertex.
            Halfedge_handle hehNewVertex = create_center_vertex(h);
            hehNewVertex->vertex()->point() = f->getRemovedVertexPos();

            // Mark all the created edges as new.
            Halfedge_around_vertex_circulator Hvc = hehNewVertex->vertex_begin();
            Halfedge_around_vertex_circulator Hvc_end = Hvc;
            CGAL_For_all(Hvc, Hvc_end)
            {
                Hvc->setNew();
                Hvc->opposite()->setNew();
            }
        }
    }
}


/**
  * Remove all the marked edges.
  */
void HiMesh::removeInsertedEdges()
{
    for (HiMesh::Halfedge_iterator hit = halfedges_begin();
         hit!=halfedges_end(); ++hit)
    {
        if(hit->isAdded()){
        	join_facet(hit);
        }
    }
}

}
