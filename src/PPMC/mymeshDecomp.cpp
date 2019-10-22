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

#include "../PPMC/configuration.h"
#include "../PPMC/frenetRotation.h"
#include "../PPMC/mymesh.h"


/**
  * Start the next decompression operation.
  */
void MyMesh::startNextDecompresssionOp()
{
    if ((float)i_curOperationId / (i_nbQuantizations + i_nbDecimations) * 100 >= i_decompPercentage)
    {

        for (MyMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
            hit->resetState();

        for (MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
            fit->resetState();

        operation = Idle;
        b_jobCompleted = true;
    } else {
        // Start the decoder.
        start_decoding(&rangeCoder);

        beginUndecimationConquest();
    }
}


/**
  * Begin an undecimation conquest.
  */
void MyMesh::beginUndecimationConquest()
{

    for (MyMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
        hit->resetState();

    for (MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
        fit->resetState();

    // Add the first halfedge to the queue.
    pushHehInit();

    // Read the min values and the ranges.
    uint16_t i16_min;
    i16_min = decode_short(&rangeCoder);
    alphaBetaMin = *(int16_t *)&i16_min;
    unsigned alphaBetaRange = decode_short(&rangeCoder);

    // Init the range coder models.
    initqsmodel(&alphaBetaModel, alphaBetaRange, 18, 1 << 17, NULL, 0);
    initqsmodel(&connectModel, 2, 10, 1 << 9, NULL, 0);

    // Set the current operation.
    operation = UndecimationConquest;
}


/**
  * One undecimation step.
  */
void MyMesh::undecimationStep()
{
    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        Face_handle f = h->facet();

        // If the face is already processed, pick the next halfedge:
        if (f->isConquered())
            continue;

        // Decode the face symbol.
        int syfreq, ltfreq;
        ltfreq = decode_culshift(&rangeCoder, 10);
        unsigned sym = qsgetsym(&connectModel, ltfreq);
        qsgetfreq(&connectModel, sym, &syfreq, &ltfreq);
        decode_update(&rangeCoder, syfreq, ltfreq, 1 << 10);
        // Update the model.
        qsupdate(&connectModel, sym);

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h;
        do
        {
            Halfedge_handle hOpp = hIt->opposite();
            assert(!hOpp->is_border());
            if (!hOpp->facet()->isConquered())
                gateQueue.push(hOpp);
            hIt = hIt->next();
        }
        while (hIt != h);

        // Decode the geometry symbol.
        if (sym == 1)
            decodeGeometrySym(h, f);
        else
            f->setUnsplittable();

        return;
    }

    // Stop the decoder.
    done_decoding(&rangeCoder);

    // Delete the models.
    deleteqsmodel(&connectModel);
    deleteqsmodel(&alphaBetaModel);

    beginInsertedEdgeDecoding();
}


/**
  * Begin the inserted edge decoding conquest.
  */
void MyMesh::beginInsertedEdgeDecoding()
{
    // Add the first halfedge to the queue.
    pushHehInit();
    operation = InsertedEdgeDecoding;

    // Init the range coder models.
    initqsmodel(&connectModel, 2, 10, 1 << 9, NULL, 0);

    // Start the decoder.
    start_decoding(&rangeCoder);
}


/**
  * One step of the inserted edge coding conquest.
  */
void MyMesh::InsertedEdgeDecodingStep()
{
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
            if (!hIt->isProcessed() && !hIt->isNew())
                gateQueue.push(hIt);
            hIt = hIt->opposite()->next();
        }

        assert(!hIt->isNew());

        float f_edgeLen = edgeLen(h);

        // Test if there is a symbol for this edge.
        // There is no symbol if the two faces of an edge are unsplitable.
        if (h->facet()->isSplittable()
            || h->opposite()->facet()->isSplittable())
        {
            // Decode the edge symbol.
            int syfreq, ltfreq;
            ltfreq = decode_culshift(&rangeCoder, 10);
            unsigned sym = qsgetsym(&connectModel, ltfreq);
            qsgetfreq(&connectModel, sym, &syfreq, &ltfreq);
            decode_update(&rangeCoder, syfreq, ltfreq, 1 << 10);
            // Update the model.
            qsupdate(&connectModel, sym);

            // Determine if the edge is original or not.
            // Mark the edge to be removed.
            if (sym != 0)
                h->setAdded();
        }

        return;
    }

    // Stop the decoder.
    done_decoding(&rangeCoder);

    // Delete the models.
    deleteqsmodel(&connectModel);

    insertRemovedVertices();
    removeInsertedEdges();

    i_curDecimationId++; // Increment the current decimation operation id.
    i_curOperationId++;
    operation = Idle;
}


/**
  * Insert center vertices.
  */
void MyMesh ::insertRemovedVertices()
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
            Point p = getPos(getQuantizedPos(barycenter(h)) + f->getResidual());
            Halfedge_handle hehNewVertex = create_center_vertex(h);
            hehNewVertex->vertex()->point() = p;

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
void MyMesh::removeInsertedEdges()
{
    for (MyMesh::Halfedge_iterator hit = halfedges_begin();
         hit!=halfedges_end(); ++hit)
    {
        if(hit->isAdded())
            join_facet(hit);
    }
}


/**
  * Decode the geometry symbols.
  */
void MyMesh::decodeGeometrySym(Halfedge_handle heh_gate, Face_handle fh)
{

    int coord[3];
    for (unsigned i = 0; i < 3; ++i)
    {
        int syfreq, ltfreq;
		// Decode the alpha and beta symbols.
		ltfreq = decode_culshift(&rangeCoder, 18);
		unsigned sym = qsgetsym(&alphaBetaModel, ltfreq);
		qsgetfreq(&alphaBetaModel, sym, &syfreq, &ltfreq);
		decode_update(&rangeCoder, syfreq, ltfreq, 1 << 18);
		// Update the alpha and beta model.
		qsupdate(&alphaBetaModel, sym);
		// Store the value.
		coord[i] = alphaBetaMin + sym;
    }

    VectorInt correction(coord[0], coord[1], coord[2]);

    fh->setSplittable();
    fh->setResidual(correction);
}
