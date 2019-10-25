prefix=/usr/local

CPPFLAGS	= -std=c++0x -g -O2 -Wno-unused-result -fno-stack-protector
INCFLAGS	= -I ./ 
LIBS		= -lgeos -lspatialindex -lboost_program_options -lstdc++ -lgmp -lmpfr -lpthread
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


#compile all the cpp files
SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

%.o: %.cpp
	$(CXX) $(INCFLAGS) $(CPPFLAGS) -c $? -o $@

#compile additional range coder .c files
RC_SRCS := $(wildcard PPMC/rangeCoder/*.c)
RC_OBJS := $(patsubst %.c,%.o,$(RC_SRCS))

%.o: %.c
	$(CC) $(INCFLAGS) -c $? -o $@

PPMC_SRCS := $(wildcard PPMC/*.cpp)
PPMC_OBJS := $(patsubst %.cpp,%.o,$(PPMC_SRCS))

TRIANGLE_SRCS := $(wildcard triangle/*.cpp)
TRIANGLE_OBJS := $(patsubst %.cpp,%.o,$(TRIANGLE_SRCS))


all: compress decompress triangle test

compress: $(RC_OBJS) $(PPMC_OBJS) test/compress.o
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
decompress: $(RC_OBJS) $(PPMC_OBJS) test/decompress.o
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@


triangle: $(TRIANGLE_OBJS) $(PPMC_OBJS) $(RC_OBJS) test/Triangle_Compute.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) $(LIBS) $(CGALFLAGS) -o ../build/$@
	
test: test/test.o
	$(CXX) $^ $(INCFLAGS) $(CPPFLAGS) -o ../build/$@
	
.PHONY: all clean

clean:
	rm -rf ppmc *.o */*.o */*/*.o

