#ifndef TDBASE_FARM2D_SRC_HRASTER_H_
#define TDBASE_FARM2D_SRC_HRASTER_H_

#include "MyRaster.h"


namespace tdbase::farm2d {

class Hraster : public MyRaster{
	double *areas = nullptr;
	bool owns_status = false;
	bool owns_areas = false;
	bool owns_mbr = false;
public:
	Hraster() = default;
	~Hraster();
	void init(double _step_x, double _step_y, int& _dimx, int& _dimy, box *mbr, bool last_layer);
	void attach_base_storage(uint8_t *_status, double *_areas);
	double get_area(int id) const { return areas[id]; }
	double *get_areas() const { return areas; }
	void set_area(int id, double area) { areas[id] = area; }
	void print();
};


} // namespace tdbase::farm2d

#endif /* TDBASE_FARM2D_SRC_HRASTER_H_ */
