#ifndef RTREE_BUILDER_3D_HPP
#define RTREE_BUILDER_3D_HPP

#include <iostream>
#include <spatialindex/SpatialIndex.h>

/* Obtain R-tree MBR of a given geometry (for usage in R-tree indexing)  */

struct mbb_3d {
	double low[3];
	double high[3];
};
class MyVisitor : public SpatialIndex::IVisitor
{
	public:
		std::vector<SpatialIndex::id_type> matches; // contains ids of matching objects

	public:	
		MyVisitor() {}

		~MyVisitor() {
			matches.clear();
		}

		void visitNode(const SpatialIndex::INode& n) {}
		void visitData(std::string &s) {}

		void visitData(const SpatialIndex::IData& d)
		{
			matches.push_back(d.getIdentifier());
		}

		void visitData(std::vector<const SpatialIndex::IData*>& v) {}
		void visitData(std::vector<uint32_t>& v){}
};

/* Customized data stream to read from a mapping of int to geomery */
class CustomDataStream : public SpatialIndex::IDataStream
{
	public:
		CustomDataStream(std::vector<struct mbb_3d *> *inputdata ) : m_pNext(0), len(0), m_id(0)
	{
		if (inputdata->empty())
			throw Tools::IllegalArgumentException("Input size is ZERO.");
		shapes = inputdata;
		len = inputdata->size();
		iter = shapes->begin();
		readNextEntry();
	}

		virtual ~CustomDataStream()
		{
			if (m_pNext != 0) delete m_pNext;
		}

		virtual SpatialIndex::IData* getNext()
		{
			if (m_pNext == 0) return 0;

			SpatialIndex::RTree::Data* ret = m_pNext;
			m_pNext = 0;
			readNextEntry();
			return ret;
		}

		virtual bool hasNext()
		{
			return (m_pNext != 0);
		}

		virtual uint32_t size()
		{
			return len;
			//throw Tools::NotSupportedException("Operation not supported.");
		}

		virtual void rewind()
		{
			if (m_pNext != 0)
			{
				delete m_pNext;
				m_pNext = 0;
			}

			m_id  = 0;
			iter = shapes->begin();
			readNextEntry();
		}

		void readNextEntry()
		{
			if (iter != shapes->end())
			{
				//std::cerr<< "readNextEntry m_id == " << m_id << std::endl;
				SpatialIndex::Region r((*iter)->low, (*iter)->high, 3);
				m_pNext = new SpatialIndex::RTree::Data(sizeof(double), reinterpret_cast<byte*>((*iter)->low), r, m_id);
				iter++;
				m_id++;
			}
		}

		SpatialIndex::RTree::Data* m_pNext;
		std::vector<struct mbb_3d*> * shapes;
		std::vector<struct mbb_3d*>::iterator iter;

		int len;
		SpatialIndex::id_type m_id;
};

class GEOSDataStreamFileTile : public SpatialIndex::IDataStream
{
	public:
		std::map<SpatialIndex::id_type, std::string> *id_tiles;
	
	public:
		GEOSDataStreamFileTile(char *input_file, 
			std::map<SpatialIndex::id_type, std::string> *id_tiles_ptr) : m_pNext(0)
	{
		m_fin.open(input_file);
		id_tiles = id_tiles_ptr;
		m_id = 0;
		if (! m_fin)
			throw Tools::IllegalArgumentException("Input file not found.");

		readNextEntry();
	}

		virtual ~GEOSDataStreamFileTile()
		{
			if (m_pNext != 0) delete m_pNext;
		}

		virtual SpatialIndex::IData* getNext()
		{
			if (m_pNext == 0) return 0;

			SpatialIndex::RTree::Data* ret = m_pNext;
			m_pNext = 0;
			readNextEntry();
			return ret;
		}

		virtual bool hasNext()
		{
			return (m_pNext != 0);
		}

		virtual uint32_t size()
		{
			throw Tools::NotSupportedException("Operation not supported.");
		}

		virtual void rewind()
		{
			if (m_pNext != 0)
			{
				delete m_pNext;
				m_pNext = 0;
			}

			m_fin.seekg(0, std::ios::beg);
			readNextEntry();
			m_id = 0;
		}

		void readNextEntry()
		{
			std::string tile_id;
			double low[3], high[3];


			m_fin >> tile_id >> low[0] >> low[1] >> low[2]  >>  high[0] >> high[1] >> high[2];
			/* store tile_id */

			if (m_fin.good())
			{
				SpatialIndex::Region r(low, high, 3);
				m_pNext = new SpatialIndex::RTree::Data(sizeof(double), reinterpret_cast<byte*>(low), r, m_id);
				/* Use spatialproc struct */
				//stop.id_tiles[m_id] = tile_id;
				(*id_tiles)[m_id] = tile_id;
			}
			m_id++;
		}

		SpatialIndex::id_type m_id;
		std::ifstream m_fin;
		SpatialIndex::RTree::Data* m_pNext;
};

#endif

