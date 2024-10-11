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
#include <fstream>
#include <unistd.h>
#include "himesh.h"

namespace tdbase{

// Write a given number of bits in a buffer.
void writeBits(uint32_t data, unsigned i_nbBits, char *p_dest,
               unsigned &i_bitOffset, size_t &offset)
{
    assert(i_nbBits <= 25);

    uint32_t dataToAdd = data << (32 - i_nbBits - i_bitOffset);
    // Swap the integer bytes because the x86 architecture is little endian.
    dataToAdd = __builtin_bswap32(dataToAdd); // Call a GCC builtin function.

    // Write the data.
    *(uint32_t *)p_dest |= dataToAdd;

    // Update the size and offset.
    offset += (i_bitOffset + i_nbBits) / 8;
    i_bitOffset = (i_bitOffset + i_nbBits) % 8;
}


/**
  * Read a given number of bits in a buffer.
  */
uint32_t readBits(unsigned i_nbBits, char *p_src,
                  unsigned &i_bitOffset, size_t &offset)
{
    assert(i_nbBits <= 25);

    // Build the mask.
    uint32_t mask = 0;
    for (unsigned i = 0; i < 32 - i_bitOffset; ++i)
        mask |= 1 << i;
    // Swap the mask bytes because the x86 architecture is little endian.
    mask = __builtin_bswap32(mask); // Call a GCC builtin function.

    uint32_t data = *(uint32_t *)p_src & mask;

    // Swap the integer bytes because the x86 architecture is little endian.
    data = __builtin_bswap32(data); // Call a GCC builtin function.

    data >>= 32 - i_nbBits - i_bitOffset;

    // Update the size and offset.
    offset += (i_bitOffset + i_nbBits) / 8;
    i_bitOffset = (i_bitOffset + i_nbBits) % 8;

    return data;
}


// Write a floating point number in the data buffer.
void HiMesh::writeFloat(float f)
{
    *(float *)(p_data + dataOffset) = f;
    dataOffset += sizeof(float);
}


/**
  * Read a floating point number in the data buffer.
  */
float HiMesh::readFloat()
{
    float f = *(float *)(p_data + dataOffset);
    dataOffset += sizeof(float);
    return f;
}

// Write a floating point number in the data buffer.
void HiMesh::writePoint(Point &p)
{
	for (unsigned i = 0; i < 3; ++i){
		writeFloat((float)p[i]);
	}
}

// Write a floating point number in the data buffer.
Point HiMesh::readPoint()
{
	float coord[3];
	for (unsigned i = 0; i < 3; ++i){
		coord[i] = readFloat();
	}
	Point pt(coord[0], coord[1], coord[2]);
	return pt;
}

/**
  * Read an integer in the data buffer.
  */
int HiMesh::readInt()
{
    int i = *(int *)(p_data + dataOffset);
    dataOffset += sizeof(int);
    return i;
}

// Write an integer in the data buffer
void HiMesh::writeInt(int i)
{
    *(int *)(p_data + dataOffset) = i;
    dataOffset += sizeof(int);
}

/**
  * Read a 16 bits integer in the data buffer.
  */
int16_t HiMesh::readInt16()
{
    int16_t i = *(int16_t *)(p_data + dataOffset);
    dataOffset += sizeof(int16_t);
    return i;
}


// Write a 16 bits integer in the data buffer
void HiMesh::writeInt16(int16_t i)
{
    *(int16_t *)(p_data + dataOffset) = i;
    dataOffset += sizeof(int16_t);
}

/**
  * Read a 16 bits integer in the data buffer.
  */
uint16_t HiMesh::readuInt16()
{
    uint16_t i = *(uint16_t *)(p_data + dataOffset);
    dataOffset += sizeof(uint16_t);
    return i;
}


// Write a 16 bits integer in the data buffer
void HiMesh::writeuInt16(uint16_t i)
{
    *(uint16_t *)(p_data + dataOffset) = i;
    dataOffset += sizeof(uint16_t);
}

/**
  * Read a byte in the data buffer.
  */
unsigned char HiMesh::readChar()
{
	unsigned char  i = *(unsigned char  *)(p_data + dataOffset);
    dataOffset += sizeof(unsigned char );
    return i;
}

// Write a byte in the data buffer
void HiMesh::writeChar(unsigned char  i)
{
    *(unsigned char *)(p_data + dataOffset) = i;
    dataOffset += sizeof(unsigned char );
}

}
