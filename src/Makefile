prefix=/usr/local
CC = gcc
CXX = g++
NVCC = nvcc

CPPFLAGS	= -std=c++14 -g -O2 -Wno-unused-result -DBOOST_ALLOW_DEPRECATED_HEADERS -DBOOST_BIND_GLOBAL_PLACEHOLDERS -DCGAL_HAS_THREADS -DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON -fopenmp -frounding-math
NVCCFLAGS	= -std=c++14 -g -O2
INCFLAGS	= -I./ -I./include -I/usr/include -I/usr/local/include
LIBS		= -L/usr/lib/x86_64-linux-gnu/ -lgmp -lmpfr -lpthread -lboost_program_options -lstdc++ -lz

CGALFLAGS	= -Wl,-rpath -frounding-math 

ifndef DEBUG
    CPPFLAGS += -DNDEBUG
endif

ifdef USE_GPU
	CPPFLAGS += -DUSE_GPU
	LIBS  += -L/usr/local/cuda/lib64 -lcuda -lcudart
endif

ifdef CGAL_MCFSKEL_DEBUG
	CPPFLAGS += -DCGAL_MCFSKEL_DEBUG -DCGAL_MCFSKEL_INFO
endif 

ifeq ($(OS),Windows_NT)
	USE_GPU = TRUE
	INCFLAGS	= -I.\include -I C:\Boost\include\boost-1_84
	INCFLAGS	+= -I C:\Users\Lenovo\Documents\GitHub\vcpkg\packages\cgal_x64-windows\include 
	INCFLAGS	+= -I C:\Users\Lenovo\Documents\GitHub\vcpkg\packages\eigen3_x64-windows\include 
	INCFLAGS	+= -I C:\Users\Lenovo\Documents\GitHub\vcpkg\packages\zlib_x64-windows\include
	INCFLAGS	+= -I "$(CUDA_PATH)\include"
	LIBS		= -L C:\Users\Lenovo\Documents\GitHub\vcpkg\packages\zlib_x64-windows\lib 
	LIBS		+= -L C:\Users\Lenovo\Documents\GitHub\vcpkg\packages\gmp_x64-windows\lib 
	LIBS		+= -L C:\Users\Lenovo\Documents\GitHub\vcpkg\packages\mpfr_x64-windows\lib 
	LIBS		+= -L C:\Boost\lib
	LIBS		+= -lgmp -lmpfr -lpthread -lstdc++ -lzlib# -lboost_program_options   
	CXX = C:\TDM-GCC-64\bin\g++
	CPPFLAGS += -DUSE_GPU
	LIBS  += -L "$(CUDA_PATH)\lib\x64" -lcuda -lcudart
else
	CPPFLAGS += -DBOOST_HAS_THREADS
endif

#compile all the cpp files
SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

%.o: %.cpp
	$(CXX) $(INCFLAGS) $(CPPFLAGS) -c $? -o $@

#compile all the cu files
SRCS_CU := $(wildcard geometry/*.cu)
OBJS_CU := $(patsubst %.cu,%_cu.o,$(SRCS_CU))
%_cu.o: %.cu
	$(NVCC) $(INCFLAGS) $(NVCCFLAGS) -c $? -o $@
		
ifndef USE_GPU
OBJS_CU = 
endif

%.o: %.c
	$(CC) $(INCFLAGS) -c $? -o $@

#SPATIAL source and object
SPATIAL_SRCS := $(wildcard spatial/*.cpp)
SPATIAL_OBJS := $(patsubst %.cpp,%.o,$(SPATIAL_SRCS))

#SPATIAL source and object
GEOMETRY_SRCS := $(wildcard geometry/*.cpp)
GEOMETRY_OBJS := $(patsubst %.cpp,%.o,$(GEOMETRY_SRCS))

GEOMETRY_OBJS += $(OBJS_CU)

#INDEX source and object
INDEX_SRCS := $(wildcard index/*.cpp)
INDEX_OBJS := $(patsubst %.cpp,%.o,$(INDEX_SRCS))

#STORAGE source and object
STORAGE_SRCS := $(wildcard storage/*.cpp)
STORAGE_OBJS := $(patsubst %.cpp,%.o,$(STORAGE_SRCS))

#JOIN source and object
JOIN_SRCS := $(wildcard join/*.cpp)
JOIN_OBJS := $(patsubst %.cpp,%.o,$(JOIN_SRCS))

#MYMESH source and object
MYMESH_SRCS := $(wildcard mymesh/*.cpp)
MYMESH_OBJS := $(patsubst %.cpp,%.o,$(MYMESH_SRCS))

all:  join generator 3dpro
	
join: tools/joiner.o $(GEOMETRY_OBJS) $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
generator: tools/generator.o $(SPATIAL_OBJS) $(INDEX_OBJS) $(GEOMETRY_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(CPPFLAGS) $(LIBS) -o ../build/$@
	
nuclei_generator: tools/generate_nuclei.o $(SPATIAL_OBJS) $(INDEX_OBJS) $(GEOMETRY_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(CPPFLAGS) $(LIBS) -o ../build/$@

3dpro: tools/3dpro_toolkit.o $(GEOMETRY_OBJS) $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

develop: tools/develop.o $(MYMESH_OBJS)
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -o ../build/$@
	
h5tooff: tools/H5ToOFF.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -lhdf5 -o ../build/$@

wrap: tools/wrap.o
	$(CXX) $^ $(CPPFLAGS) $(LIBS) -o ../build/$@
		
.PHONY: all clean

clean:
ifeq ($(OS),Windows_NT)
	del /S *.o
else
	rm -rf *.o */*.o */*/*.o
endif

