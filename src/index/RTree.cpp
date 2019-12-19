/*
 * RTree.cpp
 *
 *  Created on: Oct 30, 2019
 *      Author: teng
 */

#include "index.h"


namespace hispeed{

bool build_index_geoms(std::vector<aab_d> &geom_mbbs, SpatialIndex::ISpatialIndex* &spidx,
		SpatialIndex::IStorageManager* &storage) {
	// build spatial index on tile boundaries
	SpatialIndex::id_type  indexIdentifier;
	CustomDataStream stream(&geom_mbbs);
	storage = SpatialIndex::StorageManager::createNewMemoryStorageManager();
	spidx   = SpatialIndex::RTree::createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, stream, *storage,
			FillFactor,
			IndexCapacity,
			LeafCapacity,
			3,
			SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
	// Error checking
	return spidx->isIndexValid();
}

CustomDataStream::CustomDataStream(std::vector<aab_d> *inputdata ) : m_pNext(0), len(0), m_id(0){
	if (inputdata->empty())
		throw Tools::IllegalArgumentException("Input size is ZERO.");
	shapes = inputdata;
	len = inputdata->size();
	iter = shapes->begin();
	readNextEntry();
}

void CustomDataStream::readNextEntry(){
	if (iter != shapes->end()){
		SpatialIndex::Region r(iter->low, iter->high, 3);
		m_pNext = new SpatialIndex::RTree::Data(sizeof(double),
				reinterpret_cast<byte*>(iter->low), r, m_id);
		iter++;
		m_id++;
	}
}

}





