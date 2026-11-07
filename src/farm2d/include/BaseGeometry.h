#ifndef TDBASE_FARM2D_BASEGEOMETRY_H
#define TDBASE_FARM2D_BASEGEOMETRY_H

#include "Box.h"


namespace tdbase::farm2d {

class BaseGeometry{
public:
    box* mbr = nullptr;
};


} // namespace tdbase::farm2d

#endif // TDBASE_FARM2D_BASEGEOMETRY_H