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

namespace tdbase{

/**
  * Start the next decompression operation.
  */

void HiMesh::startNextDecompresssionOp() {
    if(global_ctx.verbose>=2 && ((float)i_curDecimationId / i_nbDecimations * 100 < i_decompPercentage||i_curDecimationId == i_nbDecimations)){
    	log("decode %d:\t%.2f\%\t[%.2f, %.2f]", i_curDecimationId, (float)i_curDecimationId / i_nbDecimations * 100, getHausdorffDistance(), getProxyHausdorffDistance());
    }

    // check if the target LOD is reached
	if (i_curDecimationId * 100.0 / i_nbDecimations >= i_decompPercentage){
		if(i_curDecimationId == i_nbDecimations){
			// reset all the hausdorff distance to 0 for the highest LOD
			// as we do not have another round of decoding to set them
			for (HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
				fit->setHausdorff(0.0);
				fit->setProxyHausdorff(0.0);
			}
		}
		b_jobCompleted = true;
		return;
	}

	// 1. reset the states. note that the states of the vertices need not to be reset
    for (HiMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
        hit->resetState();
    for (HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
        fit->resetState();

	i_curDecimationId++; // increment the current decimation operation id.
    // 2. decoding the removed vertices and add to target facets
    decodeRemovedVertices();
    // 3. decoding the inserted edge and marking the ones added
    decodeInsertedEdges();
	// 4. truly insert the removed vertices
	insertRemovedVertices();
	// 5. truly remove the added edges
	removeInsertedEdges();
	// 6. decode the Hausdorff distances for all the facets in this LOD
	HausdorffDecodingStep();

}

void HiMesh::testIteration() {
    struct timeval start = get_cur_time();
    int num = 0;
    // Add the first halfedge to the queue.
    pushHehInit();
    while (!gateQueue.empty()) {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        Face_handle f = h->facet();
        num++;
        // If the face is already processed, pick the next halfedge:
        if (f->isConquered())
            continue;

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h;
        do {
            Halfedge_handle hOpp = hIt->opposite();
            assert(!hOpp->is_border());
            if (!hOpp->facet()->isConquered())
                gateQueue.push(hOpp);
            hIt = hIt->next();
        } while (hIt != h);
    }

    logt("%d", start, num);
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
void HiMesh::decodeBaseMesh() {
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
    for (unsigned i = 0; i < i_nbVerticesBaseMesh; ++i) {
        Point pos = readPoint();
        p_pointDeque->push_back(pos);
    }

    // Read the face vertex indices.
    for (unsigned i = 0; i < i_nbFacesBaseMesh; ++i) {
    	int nv = readInt();
        uint32_t *f = new uint32_t[nv + 1];
        // Write in the first cell of the array the face degree.
        f[0] = nv;
        for (unsigned j = 1; j < nv + 1; ++j){
        	f[j] = readInt();
        }
        p_faceDeque->push_back(f);
    }

    // Let the builder do its job.
    //global_ctx.lock();
    MyMeshBaseBuilder<HalfedgeDS> builder(p_pointDeque, p_faceDeque);
    delegate(builder);
    //global_ctx.unlock();


    // Free the memory.
    for (unsigned i = 0; i < p_faceDeque->size(); ++i){
    	delete[] p_faceDeque->at(i);
    }
    delete p_faceDeque;
    delete p_pointDeque;

    // Read the bidirectional hausdorff distance
    for(unsigned i=0;i<i_nbDecimations;i++){
    	float proxyhausdorff = readFloat();
    	float hausdorff = readFloat();
    	globalHausdorfDistance.push_back(pair<float, float>(proxyhausdorff, hausdorff));
    	//log("decode hausdorff: %d %f %f",i, proxyhausdorff, hausdorff);
    }

    // load the Hausdorff distances for the base LOD
	HausdorffDecodingStep();
}

void HiMesh::decodeRemovedVertices(){

    // Add the first halfedge to the queue.
	pushHehInit();
	while (!gateQueue.empty()) {
		Halfedge_handle h = gateQueue.front();
		gateQueue.pop();

		Face_handle f = h->facet();

		// If the face is already processed, pick the next halfedge:
		if (f->isConquered())
			continue;
		// Add the other halfedges to the queue
		Halfedge_handle hIt = h;
		do {
			Halfedge_handle hOpp = hIt->opposite();
			assert(!hOpp->is_border());
			if (!hOpp->facet()->isConquered())
				gateQueue.push(hOpp);
			hIt = hIt->next();
		} while (hIt != h);

		// Decode the face symbol.
		unsigned sym = readChar();
		if (sym == 1){
			Point rmved = readPoint();
			f->setSplittable();
			f->setRemovedVertexPos(rmved);
		} else {
			f->setUnsplittable();
		}
	}
}

/**
  * One step of the inserted edge coding conquest.
  */
void HiMesh::decodeInsertedEdges() {
    pushHehInit();
    while (!gateQueue.empty()) {
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
        if (h->facet()->isSplittable() || h->opposite()->facet()->isSplittable()) {
            // Decode the edge symbol.
            unsigned sym = readChar();
            // Determine if the edge is original or not.
            // Mark the edge to be removed.
            if (sym != 0)
                h->setAdded();
        }

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h->next();
        while (hIt->opposite() != h) {
            if (!hIt->isProcessed() && !hIt->isNew())
                gateQueue.push(hIt);
            hIt = hIt->opposite()->next();
        }
        assert(!hIt->isNew());
    }
}

void HiMesh::HausdorffDecodingStep(){

	if(i_curDecimationId >= i_nbDecimations){
		return;
	}
	//log("DecimationId: %d", this->i_curDecimationId);
	int idx = 0;
	for(HiMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
		fit->resetHausdorff();
		unsigned char hausdorff_code = readChar();
		unsigned char proxyhausdorff_code = readChar();
		float hausdorff = readFloat();
		float proxyhausdorff = readFloat();

		if(use_byte_coding){
			// decode the hausdorf distance symbols
			//printf("%f %f -> ", hausdorff, proxyhausdorff);
			hausdorff = hausdorff_code * getHausdorffDistance()/127.0;
			proxyhausdorff = proxyhausdorff_code * getProxyHausdorffDistance()/127.0;
			//printf("%f %f (%f %f)\n", hausdorff, proxyhausdorff, getHausdorffDistance(), getProxyHausdorffDistance());
		}

		fit->setHausdorff(hausdorff);
		fit->setProxyHausdorff(proxyhausdorff);
		if(global_ctx.verbose>=3)
		{
			log("decode face %d:\t%.2f %.2f", idx++, proxyhausdorff, hausdorff);
		}
	}
}

/**
  * Insert center vertices.
  */
pthread_mutex_t mtx;

void HiMesh::insertRemovedVertices() {

    // Add the first halfedge to the queue.
    pushHehInit();
    while (!gateQueue.empty()) {
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
        do {
            Halfedge_handle hOpp = hIt->opposite();
            assert(!hOpp->is_border());
            if (!hOpp->facet()->isProcessed())
                gateQueue.push(hOpp);
            hIt = hIt->next();
        } while (hIt != h);
        assert(!h->isNew());  
        
        if (f->isSplittable()) {
            pthread_mutex_lock(&mtx);
            // Insert the vertex.
            Halfedge_handle hehNewVertex = create_center_vertex(h);
            pthread_mutex_unlock(&mtx);

            hehNewVertex->vertex()->point() = f->getRemovedVertexPos();
            // Mark all the created edges as new.
            Halfedge_around_vertex_circulator Hvc = hehNewVertex->vertex_begin();
            Halfedge_around_vertex_circulator Hvc_end = Hvc;
            CGAL_For_all(Hvc, Hvc_end) {
                Hvc->setNew();
                Hvc->opposite()->setNew();
            }
        }
    }
}

/**
  * Remove all the marked edges.
  */
void HiMesh::removeInsertedEdges() {
    pthread_mutex_lock(&mtx);
    // todo: use locking for avoid competetion in CGAL
	for (HiMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); hit++) {
        if(hit->isAdded()){
            join_facet(hit);
		}
	}
    pthread_mutex_unlock(&mtx);
}

}
