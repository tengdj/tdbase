prefix=/usr/local

CPPFLAGS	= -std=c++14 -g -O2 -Wno-unused-result
NVCCFLAGS	= -std=c++14 -g -O2
INCFLAGS	= -I ./ -I /usr/include -I ./include
LIBS		= -L/usr/local/cgal/lib -L/usr/lib/x86_64-linux-gnu/ -lgmp -lmpfr -lpthread -lgeos -lboost_program_options -lstdc++ -lspatialindex -lz

SLIBS		= -L/usr/local/cgal/lib -L/usr/lib/x86_64-linux-gnu/ -lgmp -lmpfr -lpthread -lboost_program_options -lz -lstdc++
CGALFLAGS	= -Wl,-rpath

ifdef DEBUG
    CPPFLAGS += -DNDEBUG
else
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

#compile additional range coder .c files
RC_SRCS := $(wildcard PPMC/rangeCoder/*.c)
RC_OBJS := $(patsubst %.c,%.o,$(RC_SRCS))

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

#INDEX source and object
INDEX_SRCS := $(wildcard index/*.cpp)
INDEX_OBJS := $(patsubst %.cpp,%.o,$(INDEX_SRCS))

#STORAGE source and object
STORAGE_SRCS := $(wildcard storage/*.cpp)
STORAGE_OBJS := $(patsubst %.cpp,%.o,$(STORAGE_SRCS))

#JOIN source and object
JOIN_SRCS := $(wildcard join/*.cpp)
JOIN_OBJS := $(patsubst %.cpp,%.o,$(JOIN_SRCS))

all: compress join generator fixer h5tooff 3dpro

compress: tools/compress.o $(SPATIAL_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
join: tools/joiner.o $(GEOMETRY_OBJS) $(OBJS_CU) $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS) $(RC_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
generator: tools/data_generator.o $(SPATIAL_OBJS) $(INDEX_OBJS) $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(CPPFLAGS) $(LIBS) -o ../build/$@

fixer: tools/fixer.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -o ../build/$@
	
h5tooff: tools/H5ToOFF.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -lhdf5 -o ../build/$@

3dpro: tools/3dpro_toolkit.o $(GEOMETRY_OBJS) $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS) $(RC_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
cgal: tools/cgal.o $(GEOMETRY_OBJS) $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS) $(RC_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

.PHONY: all clean

clean:
	rm -rf *.o */*.o */*/*.o

