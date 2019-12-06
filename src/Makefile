prefix=/usr/local

CPPFLAGS	= -std=c++14 -g -O2 -Wno-unused-result
NVCCFLAGS	= -std=c++14 -g -O2
INCFLAGS	= -I ./ -I /usr/local/cgal/include
LIBS		= -L/usr/local/cgal/lib -L/usr/local/cuda/lib64 -L/usr/lib/x86_64-linux-gnu/ \
-lgmp -lmpfr -lpthread -lcuda -lcudart -lgeos -lboost_program_options -lstdc++ 
CGALFLAGS	= -lCGAL -Wl,-rpath

ifdef DEBUG
    CPPFLAGS += -DDEBUG
else
    CPPFLAGS += -DNDEBUG
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
SRCS_CU := $(wildcard */*.cu)
OBJS_CU := $(patsubst %.cu,%_cu.o,$(SRCS))
%_cu.o: %.cu
	$(NVCC) $(INCFLAGS) $(NVCCFLAGS) -c $? -o $@

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

#PARTITION source and object
PARTITION_SRCS := $(wildcard partition/*.cpp)
PARTITION_OBJS := $(patsubst %.cpp,%.o,$(PARTITION_SRCS))

#INDEX source and object
INDEX_SRCS := $(wildcard index/*.cpp)
INDEX_OBJS := $(patsubst %.cpp,%.o,$(INDEX_SRCS))

#STORAGE source and object
STORAGE_SRCS := $(wildcard storage/*.cpp)
STORAGE_OBJS := $(patsubst %.cpp,%.o,$(STORAGE_SRCS))

#JOIN source and object
JOIN_SRCS := $(wildcard join/*.cpp)
JOIN_OBJS := $(patsubst %.cpp,%.o,$(JOIN_SRCS))

all: compress decompress partition getoff join inter generator

compress: test/compress.o $(SPATIAL_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
decompress: test/decompress.o $(SPATIAL_OBJS) $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

getoff: test/getoff.o $(SPATIAL_OBJS) $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

partition: test/partitioner.o $(PARTITION_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS) $(RC_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
join: test/joiner.o geometry/TriInt.o geometry/SegDist.o geometry/SegDist_cu.o $(JOIN_OBJS) $(STORAGE_OBJS) $(INDEX_OBJS) $(SPATIAL_OBJS) $(RC_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

inter: test/intersection.o geometry/TriInt.o 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

generator: test/data_generator.o $(SPATIAL_OBJS) $(INDEX_OBJS) $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
test: test/test.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -o ../build/$@

	
.PHONY: all clean

clean:
	rm -rf ppmc *.o */*.o */*/*.o

