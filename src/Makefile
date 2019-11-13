prefix=/usr/local

CPPFLAGS	= -std=c++11 -g -O2
INCFLAGS	= -I ./ 
LIBS		= -lgeos -lspatialindex -lboost_program_options -lstdc++ -L/usr/local/cuda/lib64 -lgmp -lmpfr -lpthread -lcuda -lcudart
CGALFLAGS	= -lCGAL -Wl,-rpath
NVCCFLAGS 	= -arch sm_30

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
	$(NVCC) $(INCFLAGS) $(CPPFLAGS) $(NVCCFLAGS) -c $? -o $@

#compile additional range coder .c files
RC_SRCS := $(wildcard PPMC/rangeCoder/*.c)
RC_OBJS := $(patsubst %.c,%.o,$(RC_SRCS))

%.o: %.c
	$(CC) $(INCFLAGS) -c $? -o $@

PPMC_SRCS := $(wildcard PPMC/*.cpp)
PPMC_OBJS := $(patsubst %.cpp,%.o,$(PPMC_SRCS))

TRIANGLE_SRCS := $(wildcard triangle/*.cpp)
TRIANGLE_OBJS := $(patsubst %.cpp,%.o,$(TRIANGLE_SRCS))


all: compress decompress triangle test skeleton

compress: test/compress.o $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
decompress: test/decompress.o $(RC_OBJS) $(PPMC_OBJS) 
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

skeleton: test/skeleton.o index/Skeleton.o $(RC_OBJS) $(PPMC_OBJS)
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@



triangle:  test/Triangle_Compute.o triangle/TriDist_cu.o triangle/TriDist.o $(PPMC_OBJS) $(RC_OBJS)
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) $(LIBS) $(CGALFLAGS) -o ../build/$@
	
test: test/test.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -o ../build/$@
	

	
.PHONY: all clean

clean:
	rm -rf ppmc *.o */*.o */*/*.o

