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

#ifndef MYMESHBASEBUILDER_H
#define MYMESHBASEBUILDER_H

#include <CGAL/Polyhedron_incremental_builder_3.h>


template <class HDS> class MyMeshBaseBuilder : public CGAL::Modifier_base<HDS>
{
public:
    MyMeshBaseBuilder(std::deque<Point> *p_pointDeque, std::deque<uint32_t *> *p_faceDeque)
        : p_pointDeque(p_pointDeque), p_faceDeque(p_faceDeque) {}

    void operator()(HDS& hds)
    {
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        size_t nbVertices = p_pointDeque->size();
        size_t nbFaces = p_faceDeque->size();

        B.begin_surface(nbVertices, nbFaces);

        for (unsigned i = 0; i < nbVertices; ++i)
            B.add_vertex(p_pointDeque->at(i));

        for (unsigned i = 0; i < nbFaces; ++i)
        {
            B.begin_facet();
            uint32_t *f = p_faceDeque->at(i);
            for (unsigned j = 1; j < f[0] + 1; ++j)
                B.add_vertex_to_facet(f[j]);
            B.end_facet();
        }

        B.end_surface();
    }

private:
    std::deque<Point> *p_pointDeque;
    std::deque<uint32_t *> *p_faceDeque;
};


#endif // MYMESHBASEBUILDER_H
