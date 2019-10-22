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

#ifndef FRENETROTATION_H
#define FRENETROTATION_H

#include <CGAL/Cartesian/MatrixC33.h>
#include "../PPMC/mymesh.h"

// Double precision needed for the frenet rotation computations.
typedef CGAL::MatrixC33<MyKernelDouble> Matrix;

void determineFrenetFrame(const MyMesh::Halfedge_handle &heh_gate,
                          const Vector &normal,
                          Vector &t1, Vector &t2);
VectorInt frenetRotation(const VectorInt &Dist, const Vector &T1,
                         const Vector &T2, const Vector &normal);
VectorInt invFrenetRotation(const VectorInt &Frenet, const Vector &T1,
                            const Vector &T2, const Vector &normal);

#endif // FRENETROTATION_H
