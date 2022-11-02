prefix=/usr/local

CPPFLAGS	= -std=c++14 -g -O2 -Wno-unused-result
NVCCFLAGS	= -std=c++14 -g -O2
INCFLAGS	= -I./ -I./include -I/usr/include -I/usr/local/include -I/usr/local/include/eigen3 -I/usr/local/cgal/include
LIBS		= -L/usr/local/cgal/lib -L/usr/lib/x86_64-linux-gnu/ -lgmp -lmpfr -lpthread -lboost_program_options -lstdc++ -lz

SLIBS		= -L/usr/local/cgal/lib -L/usr/lib/x86_64-linux-gnu/ -lgmp -lmpfr -lpthread -lboost_program_options -lz -lstdc++
CGALFLAGS	= -Wl,-rpath -frounding-math

ifndef DEBUG
    CPPFLAGS += -DNDEBUG
endif

ifdef USE_GPU
	CPPFLAGS += -DUSE_GPU
	LIBS  += -L/usr/local/cuda/lib64 -lcuda -lcudart
endif

ifdef CGAL_MCFSKEL_DEBUG
	CPPFLAGS += -DCGAL_MCFSKEL_DEBUG
endif 
 
   
CC = gcc
CXX = g++
NVCC = nvcc


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
	
#PPMC source and object
PPMC_SRCS := $(wildcard PPMC/*.cpp)
PPMC_OBJS := $(patsubst %.cpp,%.o,$(PPMC_SRCS))

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

#compress fixer h5tooff

compress: tools/compress.o $(SPATIAL_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
join: tools/joiner.o $(GEOMETRY_OBJS) $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
generator: tools/generator.o $(SPATIAL_OBJS) $(INDEX_OBJS) $(PPMC_OBJS) $(GEOMETRY_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(CPPFLAGS) $(LIBS) -o ../build/$@
	
fixer: tools/fixer.o $(MYMESH_OBJS)
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -o ../build/$@
	
h5tooff: tools/H5ToOFF.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -lhdf5 -o ../build/$@

3dpro: tools/3dpro_toolkit.o $(GEOMETRY_OBJS) $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
.PHONY: all clean

clean:
	rm -rf *.o */*.o */*/*.o

