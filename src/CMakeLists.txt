cmake_minimum_required(VERSION 3.25)
 
project(TDBase)

# set (CMAKE_TOOLCHAIN_FILE "C:/Users/Lenovo/Documents/GitHub/vcpkg/scripts/buildsystems/vcpkg.cmake")
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
 
find_package(OpenMesh)
find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(ZLIB REQUIRED)
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

include_directories(include) 

if (OpenMesh_FOUND)
	add_executable (openmesh simplifier/openmesh.cpp)
	target_link_libraries(openmesh PRIVATE OpenMeshCore OpenMeshTools gmp)
endif()

add_executable (gauland_heckbert simplifier/gauland_heckbert.cpp)
target_link_libraries(gauland_heckbert PRIVATE gmp)

add_executable (lindstrom_turk simplifier/lindstrom_turk.cpp)
target_link_libraries(lindstrom_turk PRIVATE gmp)

if(DEFINED ENV{USE_GPU})
    message (STATUS "cuda")
    enable_language("CUDA")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_GPU")
    add_executable (tdbase tools/tdbase.cpp ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${GEOMETRY_CUDA_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
    target_link_libraries(tdbase PRIVATE gmp mpfr ${Boost_LIBRARIES} ZLIB::ZLIB pthread OpenMP::OpenMP_CXX cuda)
    
    add_executable (simulator tools/simulator.cpp ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${GEOMETRY_CUDA_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
	target_link_libraries(simulator PRIVATE gmp mpfr ${Boost_LIBRARIES} ZLIB::ZLIB pthread OpenMP::OpenMP_CXX cuda)
else()
    message (STATUS "non-cuda")
    add_executable (tdbase tools/tdbase.cpp ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
    target_link_libraries(tdbase PRIVATE gmp mpfr ${Boost_LIBRARIES} ZLIB::ZLIB pthread OpenMP::OpenMP_CXX)
    add_executable (simulator tools/simulator.cpp ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
	target_link_libraries(simulator PRIVATE gmp mpfr ${Boost_LIBRARIES} ZLIB::ZLIB pthread OpenMP::OpenMP_CXX)
endif()
    
add_executable (wrap simplifier/wrap.cpp)
target_link_libraries(wrap PRIVATE gmp mpfr)

add_executable (mymesh ${MYMESH_SRCS})
