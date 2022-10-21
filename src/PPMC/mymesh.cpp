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
#include <algorithm>


MyMesh::MyMesh(unsigned i_decompPercentage,
               const int i_mode,
			   const char* data,
			   long length) :
    CGAL::Polyhedron_3< CGAL::Simple_cartesian<float>, MyItems >(), i_mode(i_mode),
    b_jobCompleted(false), i_curDecimationId(0), dataOffset(0),
    i_decompPercentage(i_decompPercentage)
{
	assert(length>0);

	srand(PPMC_RANDOM_CONSTANT);
    if (i_mode == COMPRESSION_MODE_ID) {
        // Create the compressed data buffer.
        d_capacity = 3*length;
        p_data = new char[d_capacity];
        // Fill the buffer with 0.
        for (size_t i = 0; i < d_capacity; ++i) {
           p_data[i] = 0;
        }
	    std::istringstream is;
	    is.str(data);
	    std::istringstream os;
	    os.str(data);
        is >> *this;
        if(size_of_facets()==0){
            std::cerr<<"failed to parse the OFF file into Polyhedron"<<endl;
            exit(EXIT_FAILURE);
        }

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

		// Set the vertices of the edge that is the departure of the coding and decoding conquests.
		vh_departureConquest[0] = halfedges_begin()->opposite()->vertex();
		vh_departureConquest[1] = halfedges_begin()->vertex();
    } else {
		p_data = new char[length];
		memcpy(p_data, data, length);
    	readCompressedData();
        // Set the vertices of the edge that is the departure of the coding and decoding conquests.
        vh_departureConquest[0] = vertices_begin();
        vh_departureConquest[1] = ++vertices_begin();
    }
    i_nbVerticesInit = size_of_vertices();
    i_nbFacetsInit = size_of_facets();
}


MyMesh::~MyMesh(){
	if(p_data!=NULL){
	   delete[] p_data;
	}
}

size_t MyMesh::size_of_triangles(){
	int tri_num = 0;
	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		tri_num += f->facet_degree()-2;
	}
	return tri_num;
}

/**
  * Finish completely the current operation.
  */
void MyMesh::completeOperation()
{
    while (!b_jobCompleted)
    {
        if (i_mode == COMPRESSION_MODE_ID)
            startNextCompresssionOp();
        else{
            startNextDecompresssionOp();
        }
    }
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
