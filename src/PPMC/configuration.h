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

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

// Size of the compressed data buffer.
#define BUFFER_SIZE 10 * 1024 * 1024

#define NB_BITS_FACE_DEGREE_BASE_MESH 3

#define COMPRESSION_MODE_ID 0
#define DECOMPRESSION_MODE_ID 1

#define INV_ALPHA 2
#define INV_GAMMA 2

#define PPMC_RANDOM_CONSTANT 0315

using namespace std;

#endif // CONFIGURATION_H
