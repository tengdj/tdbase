/*****************************************************************************
* Copyright (C) 2011 Adrien Maglo and Clément Courbet
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

#include "PPMC/mymesh.h"
#include "PPMC/configuration.h"
#include "PPMC/frenetRotation.h"

#include <algorithm>


MyMesh::MyMesh(unsigned i_decompPercentage,
               const int i_mode,
               unsigned i_quantBits,
               bool b_allowConcaveFaces,
			   const char* data,
			   long length
			   ) :
    CGAL::Polyhedron_3< CGAL::Simple_cartesian<float>, MyItems >(), i_mode(i_mode),
    b_jobCompleted(false), operation(Idle),
    i_curDecimationId(0), i_curQuantizationId(0), i_curOperationId(0),
    i_levelNotConvexId(0), b_allowConvexity(b_allowConcaveFaces), connectivitySize(0),
    geometrySize(0), i_quantBits(i_quantBits), dataOffset(0),
    i_decompPercentage(i_decompPercentage)
	{
	assert(length>0);
    // Create the compressed data buffer.
    p_data = new char[length];
    // Fill the buffer with 0.
    for (size_t i = 0; i < length; ++i) {
       p_data[i] = 0;
    }

    // Initialize the range coder structure.
    rangeCoder.p_data = p_data;
    rangeCoder.p_dataOffset = &dataOffset;

    if (i_mode == COMPRESSION_MODE_ID) {
	    std::istringstream is;
	    is.str(data);
        is >> *this;
		if (keep_largest_connected_components(1) != 0){
			std::cerr << "Can't compress the mesh." << std::endl;
			std::cerr << "The codec doesn't handle meshes with several connected components." << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!is_closed()){
			std::cerr << "Can't compress the mesh." << std::endl;
			std::cerr << "The codec doesn't handle meshes with borders." << std::endl;
			exit(EXIT_FAILURE);
		}

		computeBoundingBox();
		determineQuantStep();
		quantizeVertexPositions();

		// Set the vertices of the edge that is the departure of the coding and decoding conquests.
		vh_departureConquest[0] = halfedges_begin()->opposite()->vertex();
		vh_departureConquest[1] = halfedges_begin()->vertex();
    } else {
        memcpy(p_data, data, length);
    	readCompressedData();
        // Set the vertices of the edge that is the departure of the coding and decoding conquests.
        vh_departureConquest[0] = vertices_begin();
        vh_departureConquest[1] = ++vertices_begin();

        if (i_levelNotConvexId < i_nbDecimations){
            // Decompress until the first convex LOD is reached.
            while (i_curDecimationId < i_levelNotConvexId) {
                batchOperation();
            }
            batchOperation();
        }
    }
    i_nbVerticesInit = size_of_vertices();
    i_nbFacetsInit = size_of_facets();
}


MyMesh::~MyMesh(){
	if(p_data!=NULL){
	   delete[] p_data;
	}
}

/**
  * Perform one batch of steps of the current operation.
  */
void MyMesh::batchOperation()
{
    if (b_jobCompleted)
        return;
    switch (operation)
    {
    case Idle:
        if (i_mode == COMPRESSION_MODE_ID)
            startNextCompresssionOp();
        else
            startNextDecompresssionOp();
        break;
    case DecimationConquest:
        while(operation == DecimationConquest)
            decimationStep();
        break;
    case RemovedVertexCoding:
        while(operation == RemovedVertexCoding)
            RemovedVertexCodingStep();
        break;
    case InsertedEdgeCoding:
        while(operation == InsertedEdgeCoding)
            InsertedEdgeCodingStep();
        break;
    case UndecimationConquest:
        while(operation == UndecimationConquest)
            undecimationStep();
        break;
    case InsertedEdgeDecoding:
        while(operation == InsertedEdgeDecoding)
            InsertedEdgeDecodingStep();
        break;
    default:
        break;
    }
}


/**
  * Finish completely the current operation.
  */
void MyMesh::completeOperation()
{
    while (!b_jobCompleted)
        batchOperation();
}


// Compute the mesh vertex bounding box.
void MyMesh::computeBoundingBox()
{
    std::list<Point> vList;
    for(MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
        vList.push_back(vit->point());
    MyKernel::Iso_cuboid_3 bBox = CGAL::bounding_box(vList.begin(), vList.end());

    bbMin = bBox.min();
    bbMax = bBox.max();

    Vector bbVect = bbMax - bbMin;
    f_bbVolume = bbVect.x() * bbVect.y() * bbVect.z();
}


// Get the bounding box diagonal length.
float MyMesh::getBBoxDiagonal() const
{
    return sqrt(Vector(bbMin, bbMax).squared_length());
}


/**
  * Get the bounding box center point
  */
Vector MyMesh::getBBoxCenter() const
{
    return ((bbMax - CGAL::ORIGIN)
            + (bbMin - CGAL::ORIGIN)) / 2;
}


/**
  * Determine the quantization step.
  */
void MyMesh::determineQuantStep()
{
    float f_maxRange = 0;
    for (unsigned i = 0; i < 3; ++i)
    {
        float range = bbMax[i] - bbMin[i];
        if (range > f_maxRange)
            f_maxRange = range;
    }
    f_quantStep = f_maxRange / (1 << i_quantBits);
}


// Compute and store the quantized positions of the mesh vertices.
void MyMesh::quantizeVertexPositions()
{
    unsigned i_maxCoord = 1 << i_quantBits;

    // Update the positions to fit the quantization.
    for (MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
    {
        Point p = vit->point();
        PointInt quantPoint = getQuantizedPos(p);

		// Make sure the last value is in the range.
		assert(quantPoint.x() <= i_maxCoord);
		assert(quantPoint.y() <= i_maxCoord);
		assert(quantPoint.z() <= i_maxCoord);
		/* The max value is the unique that have to to be reassigned
		   because else it would get out the range. */
		quantPoint = PointInt(quantPoint.x() == i_maxCoord ? i_maxCoord - 1 : quantPoint.x(),
							  quantPoint.y() == i_maxCoord ? i_maxCoord - 1 : quantPoint.y(),
							  quantPoint.z() == i_maxCoord ? i_maxCoord - 1 : quantPoint.z());

        Point newPos = getPos(quantPoint);
        vit->point() = newPos;
    }
}


/**
  * Quantize a position
  */
PointInt MyMesh::getQuantizedPos(Point p) const
{
    return PointInt((p.x() - bbMin.x()) / (f_quantStep * (1 << i_curQuantizationId)),
                    (p.y() - bbMin.y()) / (f_quantStep * (1 << i_curQuantizationId)),
                    (p.z() - bbMin.z()) / (f_quantStep * (1 << i_curQuantizationId)));
}


/**
  * Get a position from the quantized coordinates.
  */
Point MyMesh::getPos(PointInt p) const
{
    return Point((p.x() + 0.5) * f_quantStep * (1 << i_curQuantizationId) + bbMin.x(),
                 (p.y() + 0.5) * f_quantStep * (1 << i_curQuantizationId) + bbMin.y(),
                 (p.z() + 0.5) * f_quantStep * (1 << i_curQuantizationId) + bbMin.z());
}
