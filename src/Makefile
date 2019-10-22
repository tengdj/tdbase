prefix=/usr/local

CPPFLAGS	= -std=c++0x -g -O2 -Wno-unused-result
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
OBJS := $(filter-out PPMC/ppmc.o, $(OBJS))
OBJS := $(filter-out PPMC/decomp.o, $(OBJS))

%.o: %.cpp
	$(CXX) $(INCFLAGS) $(CPPFLAGS) -c $? -o $@

#compile additional range coder .c files
RC_SRCS := $(wildcard PPMC/rangeCoder/*.c)
RC_OBJS := $(patsubst %.c,%.o,$(RC_SRCS))

%.o: %.c
	$(CC) $(INCFLAGS) -c $? -o $@


all: ppmc decomp

ppmc: $(RC_OBJS) $(OBJS) PPMC/ppmc.o
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@
	
decomp: $(RC_OBJS) $(OBJS) PPMC/decomp.o
	$(CXX) -DCGAL_USE_GMP -DCGAL_USE_MPFR -frounding-math $^ $(INCFLAGS) $(LIBS) $(CGALFLAGS) $(CPPFLAGS) -o ../build/$@

.PHONY: all clean

clean:
	rm -rf ppmc *.o */*.o */*/*.o

