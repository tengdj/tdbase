cmake_minimum_required(VERSION 3.25)
 
project(TDBase)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
 
find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
find_package(OpenMP REQUIRED)

set (CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

file(GLOB SPATIAL_SRCS
	spatial/*.cpp
)

file(GLOB GEOMETRY_SRCS
	geometry/*.cpp
)

file(GLOB GEOMETRY_CUDA_SRCS
	geometry/*.cu
)

file(GLOB INDEX_SRCS
	index/*.cpp
)

file(GLOB STORAGE_SRCS
	storage/*.cpp
)

file(GLOB JOIN_SRCS
	join/*.cpp
)

file(GLOB MYMESH_SRCS
	mymesh/*.cpp
)

file(GLOB TDBASE_SRCS
	tools/*.cpp
)

include_directories(include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_HAS_THREADS -DCGAL_HAS_THREADS -DCGAL_DONT_USE_LAZY_KERNEL")

if(DEFINED ENV{USE_GPU})
    message (STATUS "cuda")
    set(CMAKE_CUDA_ARCHITECTURES 86)
    set(CMAKE_CUDA_COMPILER /usr/bin/nvcc)
    enable_language("CUDA")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_GPU")
    add_executable (tdbase ${TDBASE_SRCS} ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${GEOMETRY_CUDA_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
    target_link_libraries(tdbase PRIVATE gmp mpfr pthread OpenMP::OpenMP_CXX cuda)
    
else()
    message (STATUS "non-cuda")
    add_executable (tdbase ${TDBASE_SRCS} ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
    target_link_libraries(tdbase PRIVATE gmp mpfr pthread OpenMP::OpenMP_CXX)
endif()
    
add_executable (wrap simplifier/wrap.cpp)
target_link_libraries(wrap PRIVATE gmp mpfr)

add_executable (mymesh ${MYMESH_SRCS})

