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

#include "math.h"
#include "../PPMC/configuration.h"
#include "../PPMC/frenetRotation.h"
#include "../PPMC/mymesh.h"


/**
  * Start the next compression operation.
  */
void MyMesh::startNextCompresssionOp()
{
    beginDecimationConquest();
}


void MyMesh::beginDecimationConquest()
{
  //printf("Begin decimation conquest n°%u.\n", i_curDecimationId);

  for(MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
        vit->resetState();

  for(MyMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
        hit->resetState();

//  int total = 0;
//  int split = 0;
//  int nsplit = 0;
//  for (MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit){
//	if(fit->isSplittable()){
//		split += count_triangle(fit);
//	}
//	if(fit->isUnsplittable()){
//		nsplit += count_triangle(fit);
//	}
//  }
//  total = split+nsplit;
//  printf("%d,%d,%d\n",total,split,nsplit);

  for(MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
        fit->resetState();


  // Select the first gate to begin the decimation.
  // teng: we always start from the middle
  // size_t i_heInitId = (float)rand() / RAND_MAX * size_of_halfedges();
  size_t i_heInitId = size_of_halfedges()/2;
  Halfedge_iterator hitInit = halfedges_begin();
  for (unsigned i = 0; i < i_heInitId; ++i)
      ++hitInit;

  hitInit->setInQueue();
  gateQueue.push((Halfedge_handle)hitInit);

  // Reset the number of removed vertices.
  i_nbRemovedVertices = 0;

  // Set the current operation.
  operation = DecimationConquest;
}


// One decimation step.
void MyMesh::decimationStep()
{

    //choose a halfedge that can be processed:
    while(!gateQueue.empty())
    {
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
        int numSplittableVerts = 0;
        Halfedge_handle unconqueredVertexHE;

        for(Halfedge_handle hh = h->next(); hh!=h; hh=hh->next())
        {
            if(isRemovable(hh->vertex()))
            {
                if(numSplittableVerts==0)
                    unconqueredVertexHE = hh;
                ++numSplittableVerts;
            }
        }

        //if all face vertices are conquered, then the current face is a null patch:
        if(numSplittableVerts==0)
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
            return;
        }
        //Is there a unique possible choice among the face vertices ?
        else //if (numSplittableVerts==1)
        {
            //in that case, cornerCut that vertex.
            h->removeFromQueue();
            vertexCut(unconqueredVertexHE);
            return;
        }
    }

    if (i_nbRemovedVertices == 0)
    {
		for(MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
		{
			if(isRemovable(vit))
				printf("Still a vertex that can be removed !\n");
		}


		operation = Idle;
		b_jobCompleted = true;
		i_curDecimationId--;
		writeCompressedData();
    }
    else
    {
        // Determine the residuals.
        determineResiduals();

        operation = RemovedVertexCoding;
        beginRemovedVertexCodingConquest();
    }

    /*
     * record the maximum volume change
     * */
    if(!b_jobCompleted){
        //log("%f", tmpMaximumcut);
    	maximumCut.push_back(tmpMaximumcut);
    }
    // reset
    tmpMaximumcut = 0;
}


/**
  * Perform the re-edging and the vertex removal.
  * Store the position of the removed vertex.
  */
MyMesh::Halfedge_handle MyMesh::vertexCut(Halfedge_handle startH)
{
        Vertex_handle v = startH->vertex();

        //make sure that the center vertex can be removed
        //assert(!v->isConquered());
        assert(v->vertex_degree()>2);

        Halfedge_handle h = startH->opposite(), end(h);
        do
        {
                assert(!h->is_border());
                Face_handle f = h->facet();
                assert(!f->isConquered()); //we cannot cut again an already cut face, or a NULL patch

                //if the face is not a triangle, cut the corner
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
                }

                //mark the vertex as conquered
                h->vertex()->setConquered();
        }
        while((h=h->opposite()->next()) != end);

        //copy the position of the center vertex:
        Point vPos = startH->vertex()->point();

        //remove the center vertex
        Halfedge_handle hNewFace = erase_center_vertex(startH);

        //now mark the new face as having a removed vertex
        hNewFace->facet()->setSplittable();
        // keep the removed vertex position.
        hNewFace->facet()->setRemovedVertexPos(vPos);

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
  * Determine the residuals to encode.
  */
void MyMesh::determineResiduals()
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

        // besides keep the residual, also update the maximum cut if needed
        if (f->isSplittable()){
        	Point rmved = f->getRemovedVertexPos();
        	Point bc = barycenter(h);
        	float cutdist = (rmved.x()-bc.x())*(rmved.x()-bc.x())+
        			(rmved.y()-bc.y())*(rmved.y()-bc.y())+
					(rmved.z()-bc.z())*(rmved.z()-bc.z());
        	tmpMaximumcut = max(tmpMaximumcut, cutdist);
        	//log("%d %f",processCount++, cutdist);
            f->setResidual(getQuantizedPos(rmved) - getQuantizedPos(bc));

            //f->setResidual(getQuantizedPos(f->getRemovedVertexPos()) - getQuantizedPos(barycenter(h)));

        }
    }
}


/**
  * Begin the removed vertex coding conquest.
  */
void MyMesh::beginRemovedVertexCodingConquest()
{
    for(MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
          fit->resetProcessedFlag();

    // Add the first halfedge to the queue.
    pushHehInit();

    // Resize the vectors to add the current conquest symbols.
    geometrySym.push_back(std::deque<VectorInt>());
    connectFaceSym.push_back(std::deque<unsigned>());

    operation = RemovedVertexCoding;
}


/**
  * One step of the removed vertex coding conquest.
  */
void MyMesh::RemovedVertexCodingStep()
{
    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        Face_handle f = h->facet();

        // If the face is already processed, pick the next halfedge:
        if (f->isProcessed())
            continue;

        // Determine face symbol.
        unsigned sym;
        bool b_split = f->isSplittable();

        // No connectivity prediction.
        sym = b_split ? 1 : 0;

        // Push the symbols.
        connectFaceSym[i_curDecimationId].push_back(sym);

        // Determine the geometry symbol.
        if (b_split)
            determineGeometrySym(h, f);

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

        return;
    }

    operation = InsertedEdgeCoding;
    beginInsertedEdgeCoding();
}


/**
  * Determine the geometry symbols.
  */
void MyMesh::determineGeometrySym(Halfedge_handle heh_gate, Face_handle fh)
{
    VectorInt distQuant = fh->getResidual();
    geometrySym[i_curDecimationId].push_back(distQuant);
}


/**
  * Begin the inserted edge coding conquest.
  */
void MyMesh::beginInsertedEdgeCoding()
{
    // Add the first halfedge to the queue.
    pushHehInit();
    // Resize the vector to add the current conquest symbols.
    connectEdgeSym.push_back(std::deque<unsigned>());

    operation = InsertedEdgeCoding;
}


/**
  * One step of the inserted edge coding conquest.
  */
void MyMesh::InsertedEdgeCodingStep()
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
            if (!hIt->isProcessed())
                gateQueue.push(hIt);
            hIt = hIt->opposite()->next();
        }

        // Don't write a symbol if the two faces of an egde are unsplitable.
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

        return;
    }

    i_curDecimationId++; // Increment the current decimation operation id.
    i_curOperationId++;
    operation = Idle;
}


/**
  * Encode an inserted edge list.
  */
void MyMesh::encodeInsertedEdges(unsigned i_operationId)
{
    // Start the encoder.
    start_encoding(&rangeCoder, 0, 0);

    std::deque<unsigned> &symbols = connectEdgeSym[i_operationId];

    assert(symbols.size() > 0);

    // Init the connectivity range coder.
    initqsmodel(&connectModel, 2, 10, 1 << 9, NULL, 1);

    unsigned i_len = symbols.size();
    for (unsigned i = 0; i < i_len; ++i)
    {
        unsigned sym = symbols[i];

        // Encode the symbol.
        int syfreq, ltfreq;
        qsgetfreq(&connectModel, sym, &syfreq, &ltfreq);
        encode_shift(&rangeCoder, syfreq, ltfreq, 10);
        qsupdate(&connectModel, sym);
    }

    unsigned i_size = done_encoding(&rangeCoder);

    connectivitySize += i_size * 8;
    // Destroy the model.
    deleteqsmodel(&connectModel);
}


/**
  * Encode the geometry and the connectivity of a removed vertex list.
  */
void MyMesh::encodeRemovedVertices(unsigned i_operationId)
{
    // Start the encoder.
    start_encoding(&rangeCoder, 0, 0);
    
    std::deque<unsigned> &connSym = connectFaceSym[i_operationId];
    std::deque<VectorInt> &geomSym = geometrySym[i_operationId];

    unsigned i_lenGeom = geomSym.size();
    unsigned i_lenConn = connSym.size();
    assert(i_lenGeom > 0);
    assert(i_lenConn > 0);

    // Determine the min and max values for the geometry coding.
    alphaBetaMin = geomSym[0].x();
    int alphaBetaMax = geomSym[0].x();

    for (unsigned i = 0; i < i_lenGeom; ++i)
    {
        VectorInt v = geomSym[i];

        if(v.x() < alphaBetaMin)
            alphaBetaMin = v.x();
        if(v.y() < alphaBetaMin)
            alphaBetaMin = v.y();
        if(v.z() < alphaBetaMin)
            alphaBetaMin = v.z();

        if(v.x() > alphaBetaMax)
            alphaBetaMax = v.x();
        if(v.y() > alphaBetaMax)
            alphaBetaMax = v.y();
        if(v.z() > alphaBetaMax)
            alphaBetaMax = v.z();
    }

    // Check that we have at least two geometry symbols.
    unsigned alphaBetaRange = std::max(alphaBetaMax - alphaBetaMin + 1, 2);

    // Write the min value and the range.
    int16_t i16_min;
    assert(alphaBetaMin >= -(1 << 14) && alphaBetaMin <= (1 << 14));
    i16_min = alphaBetaMin;
    encode_short(&rangeCoder, *(uint16_t *)&i16_min);
    assert(alphaBetaRange < (1 << 14));
    encode_short(&rangeCoder, alphaBetaRange);

    // Range coder to only measure the size of the connectivity data.
    size_t *p_dataOffsetMes = new size_t;
    *p_dataOffsetMes = 0;
    char *p_dataMes = new char[BUFFER_SIZE];
    rangecoder rangeCoderMes;
    rangeCoderMes.p_data = p_dataMes;
    rangeCoderMes.p_dataOffset = p_dataOffsetMes;
    start_encoding(&rangeCoderMes, 0, 0);

    // Init the models.
    initqsmodel(&alphaBetaModel, alphaBetaRange, 18, 1 << 17, NULL, 1);
    initqsmodel(&connectModel, 2, 10, 1 << 9, NULL, 1);

    unsigned k = 0;
    for (unsigned i = 0; i < i_lenConn; ++i)
    {
        // Encode the connectivity.
        unsigned sym = connSym[i];
        bool b_split = connSym[i];

        int syfreq, ltfreq;
        // Encode the symbol.
        qsgetfreq(&connectModel, sym, &syfreq, &ltfreq);
        encode_shift(&rangeCoder, syfreq, ltfreq, 10);
        encode_shift(&rangeCoderMes, syfreq, ltfreq, 10);
        // Update the model.
        qsupdate(&connectModel, sym);

        // Encode the geometry if necessary.
        if (b_split)
        {
            VectorInt v = geomSym[k];
            for (unsigned j = 0; j < 3; ++j)
            {

                    sym = v[j] - alphaBetaMin;
                    // Encode the alpha and beta symbols.
                    qsgetfreq(&alphaBetaModel, sym, &syfreq, &ltfreq);
                    encode_shift(&rangeCoder, syfreq, ltfreq, 18);
                    // Update the alpha and beta model.
                    qsupdate(&alphaBetaModel, sym);
            }
            k++;
        }
    }

    unsigned i_size = done_encoding(&rangeCoder);
    unsigned i_sizeConn = done_encoding(&rangeCoderMes);

    geometrySize += (i_size - i_sizeConn) * 8;
    connectivitySize += i_sizeConn * 8;

    // Destroy the models.
    deleteqsmodel(&alphaBetaModel);
    deleteqsmodel(&connectModel);

    delete[] p_dataMes;
    delete p_dataOffsetMes;
}
