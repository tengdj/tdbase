# Compatibility entry point; all commands use the integrated FARM library.
.DEFAULT_GOAL := all
.NOTPARALLEL:
BUILD_DIR ?= ../../build/cmake
GPU := $(if $(filter 1 ON on TRUE true YES yes,$(USE_GPU)),ON,OFF)
PROFILE := $(if $(filter 1 ON on TRUE true YES yes,$(PROFILE_QUERY)),ON,OFF)
TARGETS := intersect_cpu intersection_cpu within_polygon_cpu dump_polygons filter_polygons
ifeq ($(GPU),ON)
TARGETS += intersect intersection within_polygon
endif
.PHONY: all configure clean $(TARGETS)
all: configure
	cmake --build $(BUILD_DIR) --target farm2d --parallel 2
configure:
	cmake -S .. -B $(BUILD_DIR) -G Ninja -DUSE_GPU=$(GPU) -DFARM_PROFILE_QUERY=$(PROFILE) -DCMAKE_BUILD_TYPE=Release
intersect_cpu: configure
	cmake --build $(BUILD_DIR) --target farm2d_intersect_cpu --parallel 2
intersection_cpu: configure
	cmake --build $(BUILD_DIR) --target farm2d_intersection_cpu --parallel 2
within_polygon_cpu: configure
	cmake --build $(BUILD_DIR) --target farm2d_within_polygon_cpu --parallel 2
dump_polygons: configure
	cmake --build $(BUILD_DIR) --target farm2d_dump_polygons --parallel 2
filter_polygons: configure
	cmake --build $(BUILD_DIR) --target farm2d_filter_polygons --parallel 2
intersect: configure
	cmake --build $(BUILD_DIR) --target farm2d_intersect_gpu --parallel 2
intersection: configure
	cmake --build $(BUILD_DIR) --target farm2d_intersection_gpu --parallel 2
within_polygon: configure
	cmake --build $(BUILD_DIR) --target farm2d_within_polygon_gpu --parallel 2
clean:
	cmake --build $(BUILD_DIR) --target clean
