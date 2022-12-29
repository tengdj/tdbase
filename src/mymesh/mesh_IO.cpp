/*
 * mesh_parse.cpp
 *
 *  Created on: Nov 16, 2022
 *      Author: teng
 */

#include "mymesh.h"

namespace hispeed{


string Polyhedron::to_string(){

	stringstream os;
	os << "OFF" << endl;
	os << vertices.size() << " " << faces.size() << " 0\n" << endl;
	int idx = 0;
	for(Vertex *p:vertices){
		p->id = idx++;
		os << p->v[0] << " " << p->v[1] << " " << p->v[2] << endl;
	}
	for(Face *f:faces){
		os<<f->vertices.size()<<" ";
		for(Vertex *v:f->vertices){
			os<<v->id<<" ";
		}
		if(f->added){
			os<<"\t255 255 0";
		}else{
			os<<"\t0 255 0";
		}
		os<<endl;
	}
	return os.str();
}

void Polyhedron::print(){
	cout<<to_string().c_str();
}

void Polyhedron::dumpto(string fp){
	ofstream of(fp.c_str());
	of<<to_string().c_str();
	of.close();
}

//
// read vertices and faces from original OFF file
//
void Polyhedron::load(string fp) {
	string content = hispeed::read_file(fp.c_str());
	parse(content.c_str(), content.size());
}


bool parse_OFF_sign(const char *data, size_t &offset, size_t len){
	if(len<3){
		return false;
	}
	while(offset<len-3){
		if((data[offset]=='O' || data[offset]=='o')
			&&(data[offset+1]=='F' || data[offset+1]=='f')
			&&(data[offset+2]=='F' || data[offset+2]=='f')){
			offset += 3;
			return true;
		}else{
			offset++;
			continue;
		}
	}
	// no OFF is found;
	return false;
}

void skip_spaces(const char *data, size_t &offset, size_t len){
	while(offset<len){
		if(data[offset]==' '||data[offset]=='\t'||data[offset]=='\n'||data[offset]=='|'){
			offset++;
			continue;
		}
		break;
	}
}

bool parse_float(const char *data, size_t &offset, size_t len, float &ret){
	skip_spaces(data, offset, len);
	char buffer[20];
	int iter = 0;
	bool has_point = false;
	while(true){
		if(offset==len || ((data[offset]<'0' || data[offset]>'9') && data[offset]!='.')){
			if(iter>0){
				buffer[iter] = '\0';
				ret = atof(buffer);
				return true;
			}else{
				return false;
			}
		}else{
			buffer[iter++] = data[offset++];
			if(iter>=19){
				buffer[iter] = '\0';
				log("%s is too long", buffer);
				return false;
			}

			if(data[offset]=='.'){
				if(has_point){
					buffer[iter] = '\0';
					log("wrong float value", buffer);
					return false;
				}
				has_point = true;
			}
		}
	}
	return false;
}

bool parse_int(const char *data, size_t &offset, size_t len, int &ret){
	skip_spaces(data, offset, len);
	char buffer[20];
	int iter = 0;
	while(true){
		if(offset==len || data[offset]<'0' || data[offset]>'9'){
			if(iter>0){
				buffer[iter] = '\0';
				ret = atoi(buffer);
				return true;
			}else{
				return false;
			}
		}else{
			buffer[iter++] = data[offset++];
			if(iter>=19){
				buffer[iter] = '\0';
				log("%s is too long", buffer);
				return false;
			}
		}
	}
	return false;
}

bool Polyhedron::parse(string str){
	return parse(str.c_str(), str.size());
}

// parse from OFF
bool Polyhedron::parse(const char *data, size_t len){
	size_t offset = 0;
	float flt_tmp;
	int int_tmp;
	if(!parse_OFF_sign(data,offset,len)){
		log("no OFF sign found");
		return false;
	}
	int num_vertices, num_faces, num_edges;
	if(!parse_int(data, offset, len, num_vertices)||
	   !parse_int(data, offset, len, num_faces)||
	   !parse_int(data, offset, len, num_edges)){
		return false;
	}
	if(num_vertices==0){
		return true;
	}
	vector<Vertex *> vertex_array;
	vertex_array.resize(num_vertices);
	for(int i=0;i<num_vertices;i++){
		float co[3];
		for(int j=0;j<3;j++){
			if(!parse_float(data, offset, len, co[j])){
				return false;
			}
		}
		vertex_array[i] = new Vertex(co[0], co[1], co[2]);
	}
	if(num_faces==0){
		return true;
	}
	for(int i=0;i<num_faces;i++){
		int fnum = 0;
		if(!parse_int(data, offset, len, fnum) || fnum<3){
			log("invalid face numbers %d", fnum);
			return false;
		}
		vector<Vertex *> vs;
		vs.resize(fnum);
		for(int j=0;j<fnum;j++){
			int vindex = 0;
			if(!parse_int(data, offset, len, vindex) || vindex>=vertex_array.size()){
				log("invalid vertex index %d", vindex);
				return false;
			}
			vs[j] = vertex_array[vindex];
		}
		add_face(vs);
	}
	vertices.insert(vertex_array.begin(), vertex_array.end());

	return true;
}


/*
 * the operations for read/write the data
 *
 * */

/**
  * Write the compressed data to the buffer.
  */
void Polyhedron::writeCompressedData()
{

    i_nbDecimations = i_curDecimationId + 1;

    // Write the base mesh.
    writeBaseMesh();
    int i_deci = i_curDecimationId;
    assert(i_deci>0);
    while (i_deci>=0)
    {
		encodeRemovedVertices(i_deci);
		encodeInsertedEdges(i_deci);
		i_deci--;
    }
}


/**
  * Read the compressed data from the buffer.
  */
void Polyhedron::readCompressedData()
{
    // Read the base mesh.
    readBaseMesh();
}

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
void Polyhedron::writeFloat(float f)
{
    *(float *)(p_data + dataOffset) = f;
    dataOffset += sizeof(float);
}


/**
  * Read a floating point number in the data buffer.
  */
float Polyhedron::readFloat()
{
    float f = *(float *)(p_data + dataOffset);
    dataOffset += sizeof(float);
    return f;
}

/**
  * Read an integer in the data buffer.
  */
int Polyhedron::readInt()
{
    int i = *(int *)(p_data + dataOffset);
    dataOffset += sizeof(int);
    return i;
}

// Write an integer in the data buffer
void Polyhedron::writeInt(int i)
{
    *(int *)(p_data + dataOffset) = i;
    dataOffset += sizeof(int);
}

/**
  * Read a 16 bits integer in the data buffer.
  */
int16_t Polyhedron::readInt16()
{
    int16_t i = *(int16_t *)(p_data + dataOffset);
    dataOffset += sizeof(int16_t);
    return i;
}


// Write a 16 bits integer in the data buffer
void Polyhedron::writeInt16(int16_t i)
{
    *(int16_t *)(p_data + dataOffset) = i;
    dataOffset += sizeof(int16_t);
}

/**
  * Read a 16 bits integer in the data buffer.
  */
uint16_t Polyhedron::readuInt16()
{
    uint16_t i = *(uint16_t *)(p_data + dataOffset);
    dataOffset += sizeof(uint16_t);
    return i;
}


// Write a 16 bits integer in the data buffer
void Polyhedron::writeuInt16(uint16_t i)
{
    *(uint16_t *)(p_data + dataOffset) = i;
    dataOffset += sizeof(uint16_t);
}


/**
  * Read a byte in the data buffer.
  */
unsigned char Polyhedron::readChar()
{
	unsigned char  i = *(unsigned char  *)(p_data + dataOffset);
    dataOffset += sizeof(unsigned char );
    return i;
}

// Write a byte in the data buffer
void Polyhedron::writeChar(unsigned char  i)
{
    *(unsigned char *)(p_data + dataOffset) = i;
    dataOffset += sizeof(unsigned char );
}

// Write the base mesh.
void Polyhedron::writeBaseMesh()
{

    // Write the bounding box coordinate.
    for (unsigned i = 0; i < 3; ++i)
        writeFloat((float)(mbb.low[i]));
    for (unsigned i = 0; i < 3; ++i)
        writeFloat((float)(mbb.high[i]));

    // Write the number of level of decimations.
    writeInt16(i_nbDecimations);

    // Write the number of vertices and faces on 16 bits.
    unsigned i_nbVerticesBaseMesh = size_of_vertices();
    unsigned i_nbFacesBaseMesh = size_of_facets();
    writeInt(i_nbVerticesBaseMesh);
    writeInt(i_nbFacesBaseMesh);

    // Write the base mesh vertex coordinates.
    unsigned i_nbAdditionalBitsGeometry = 0;

    // Write the vertices of the edge that is the departure of the coding conquests.
    for (unsigned i = 0; i < 2; ++i)
    {
        Vertex *p = vh_departureConquest[i];
        for (unsigned j = 0; j < 3; ++j)
        {
        	writeFloat(p->v[j]);
        }
        vh_departureConquest[i]->id = i;
    }

    // Write the other vertices.
    size_t id = 2;
    for (Vertex *v:vertices)
    {
        if (v == vh_departureConquest[0] || v == vh_departureConquest[1])
            continue;
        // Write the coordinates.
        for (unsigned i = 0; i < 3; ++i)
        {
        	writeFloat(v->v[i]);
        }
        // Set an id to the vertex.
        v->id = id++;
    }

    // Write the base mesh face vertex indices.
    for (Face *fit:faces)
    {
        unsigned i_faceDegree = fit->facet_degree();
        writeInt(i_faceDegree);
        for(Vertex *v:fit->vertices){
        	writeInt(v->getId());
        }
    }

    // Write the hausdorf distance for the base faces
    for (Face *fit:faces)
    {
        writeFloat(fit->getHausdorfDistance().first);
        writeFloat(fit->getHausdorfDistance().second);
    }

    // 3dpro
    // Write the maximum volume change for each round of decimation
    assert(globalHausdorfDistance.size()==i_nbDecimations);
    for(unsigned i=0;i<i_nbDecimations;i++){
    	writeFloat(globalHausdorfDistance[i].first);
    	writeFloat(globalHausdorfDistance[i].second);
    }
}


// Read the base mesh.
void Polyhedron::readBaseMesh()
{
    // Read the bounding box min coordinate.
    for (unsigned i = 0; i < 3; ++i)
        mbb.low[i] = readFloat();

    for (unsigned i = 0; i < 3; ++i)
        mbb.high[i] = readFloat();

    // Read the number of level of detail.
    i_nbDecimations = readInt16();

    // Set the mesh bounding box.
    unsigned i_nbVerticesBaseMesh = readInt();
    unsigned i_nbFacesBaseMesh = readInt();

    vector<Vertex *> tmp_vertices;
    tmp_vertices.resize(i_nbVerticesBaseMesh);
    // Read the vertex positions.
    for (unsigned i = 0; i < i_nbVerticesBaseMesh; ++i)
    {
    	tmp_vertices[i] = new Vertex();
        for (unsigned j = 0; j < 3; ++j){
        	tmp_vertices[i]->v[j] = readFloat();
        }
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

    for (Polyhedron::Facet_iterator fit = facets_begin();
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

void Polyhedron::writeMeshOff(const char psz_filePath[])
{
    std::filebuf fb;
    fb.open(psz_filePath, std::ios::out | std::ios::trunc);
    if(fb.is_open())
    {
        std::ostream os(&fb);
        os << *this;
    }else{
    	std::cerr<<"cannot find path "<<psz_filePath<<std::endl;
    }
}


void Polyhedron::writeCurrentOperationMesh(std::string pathPrefix, unsigned i_id)
{
    // Output the current mesh in an off file.
    std::ostringstream fileName;
    fileName << pathPrefix << "_" << i_id << ".off";
    writeMeshOff(fileName.str().c_str());
}



}


