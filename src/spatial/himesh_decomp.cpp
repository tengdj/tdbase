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
    if(global_ctx.verbose>=2 && ((float)i_curDecimationId / i_nbDecimations * 100 < i_decompPercentage||i_curDecimationId == i_nbDecimations)){
    	log("decode %d:\t%.2f\%\t[%.2f, %.2f]", i_curDecimationId, (float)i_curDecimationId / i_nbDecimations * 100, getHausdorfDistance().first, getHausdorfDistance().second);
    }
	if (i_curDecimationId * 100.0 / i_nbDecimations >= i_decompPercentage){
		if(i_curDecimationId == i_nbDecimations){
			// reset all the hausdorf distance to 0 for the highest LOD
			// as we do not have another round of decoding to set them
			for (HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
				fit->setProgressive(0.0);
				fit->setConservative(0.0);
			}
		}
		b_jobCompleted = true;
		return;
	}

	// reset the states
    for (HiMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
        hit->resetState();

    for (HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
        fit->resetState();

	undecimationStep();
	InsertedEdgeDecodingStep();
	insertRemovedVertices();
	removeInsertedEdges();

	i_curDecimationId++; // increment the current decimation operation id.
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

void HiMesh::decode(int lod){
	assert(lod>=0 && lod<=100);
	assert(!this->is_compression_mode());
	if(lod<i_decompPercentage){
		return;
	}
	i_decompPercentage = lod;
	b_jobCompleted = false;
	while(!b_jobCompleted) {
		startNextDecompresssionOp();
	}
}

// Read the base mesh.
void HiMesh::readBaseMesh()
{
    // Read the bounding box
    for (unsigned i = 0; i < 3; ++i)
        mbb.low[i] = readFloat();
    for (unsigned i = 0; i < 3; ++i)
        mbb.high[i] = readFloat();

    // Read the number of level of detail.
    i_nbDecimations = readInt16();

    // Set the mesh bounding box.
    unsigned i_nbVerticesBaseMesh = readInt();
    unsigned i_nbFacesBaseMesh = readInt();

    std::deque<Point> *p_pointDeque = new std::deque<Point>();
    std::deque<uint32_t *> *p_faceDeque = new std::deque<uint32_t *>();

    // Read the vertex positions.
    for (unsigned i = 0; i < i_nbVerticesBaseMesh; ++i)
    {
        float p[3];
        for (unsigned j = 0; j < 3; ++j)
            p[j] = readFloat();
        Point pos(p[0], p[1], p[2]);
        p_pointDeque->push_back(pos);
    }

    // Read the face vertex indices.
    for (unsigned i = 0; i < i_nbFacesBaseMesh; ++i)
    {
        uint32_t *f = new uint32_t[(1 << NB_BITS_FACE_DEGREE_BASE_MESH) + 3];
        // Write in the first cell of the array the face degree.
        f[0] = readInt();
        for (unsigned j = 1; j < f[0] + 1; ++j){
        	f[j] = readInt();
        }
        p_faceDeque->push_back(f);
    }

    // Let the builder do its job.
    MyMeshBaseBuilder<HalfedgeDS> builder(p_pointDeque, p_faceDeque);
    delegate(builder);

    for (HiMesh::Facet_iterator fit = facets_begin();
         fit != facets_end(); ++fit)
    {
    	float hdist = readFloat();
    	fit->setConservative(hdist);
    	hdist = readFloat();
    	fit->setProgressive(hdist);
    }

    // Free the memory.
    for (unsigned i = 0; i < p_faceDeque->size(); ++i)
        delete[] p_faceDeque->at(i);
    delete p_faceDeque;
    delete p_pointDeque;

    // Read the maximum cutting volume
    for(unsigned i=0;i<i_nbDecimations;i++){
    	float conservative = readFloat();
    	float progressive = readFloat();
    	globalHausdorfDistance.push_back(pair<float, float>(conservative, progressive));
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
