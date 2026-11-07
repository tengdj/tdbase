cmake_minimum_required(VERSION 3.25)
 
project(TDBase)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Keep the historical USE_GPU environment variable, while allowing callers
# to explicitly choose -DUSE_GPU=ON/OFF even when that variable is exported.
set(USE_GPU_DEFAULT OFF)
if(DEFINED ENV{USE_GPU})
    set(USE_GPU_DEFAULT "$ENV{USE_GPU}")
endif()
option(USE_GPU "Build TDBase and FARM CUDA queries" ${USE_GPU_DEFAULT})

if(USE_GPU)
    if(NOT CMAKE_CUDA_COMPILER AND NOT DEFINED ENV{CUDACXX})
        find_program(TDBASE_NVCC nvcc HINTS "$ENV{CUDA_HOME}/bin" /usr/local/cuda/bin)
        if(TDBASE_NVCC)
            set(CMAKE_CUDA_COMPILER "${TDBASE_NVCC}" CACHE FILEPATH "CUDA compiler")
        endif()
    endif()
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES native)
    endif()
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
endif()
 
find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Threads REQUIRED)

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

file(GLOB TOOlS_SRCS
	tools/*.cpp
)

include_directories(include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_HAS_THREADS -DCGAL_HAS_THREADS -DCGAL_DONT_USE_LAZY_KERNEL")

if(USE_GPU)
    message (STATUS "Building TDBase and FARM with CUDA")
    add_executable (tdbase tools/tdbase.cpp tools/simulator.cpp ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${GEOMETRY_CUDA_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
    target_compile_definitions(tdbase PRIVATE USE_GPU)
    target_link_libraries(tdbase PRIVATE CUDA::cuda_driver CUDA::cudart)
else()
    message (STATUS "non-cuda")
    add_executable (tdbase tools/tdbase.cpp tools/simulator.cpp ${SPATIAL_SRCS} ${GEOMETRY_SRCS} ${INDEX_SRCS} ${STORAGE_SRCS} ${JOIN_SRCS})
endif()
    target_link_libraries(tdbase PRIVATE gmp mpfr CGAL::CGAL OpenMP::OpenMP_CXX)

add_executable(tdbase_test tools/test.cpp)

# The namespaced FARM library is linked into TDBase itself. Compatibility
# executables call the same command and API implementation.
file(GLOB FARM2D_CORE_SRCS CONFIGURE_DEPENDS farm2d/geometry/*.cpp farm2d/storage/*.cpp)
add_library(tdbase_farm2d STATIC ${FARM2D_CORE_SRCS} farm2d/engine.cpp farm2d/cli.cpp)
target_include_directories(tdbase_farm2d BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/farm2d/include)
target_include_directories(tdbase_farm2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
option(FARM_PROFILE_QUERY "Enable detailed FARM query profiling" OFF)
target_compile_definitions(tdbase_farm2d PRIVATE FARM_PROFILE_QUERY=$<BOOL:${FARM_PROFILE_QUERY}>)
target_link_libraries(tdbase_farm2d PUBLIC Boost::program_options Threads::Threads)
if(USE_GPU)
    target_sources(tdbase_farm2d PRIVATE
        farm2d/cuda/preprocessing.cu farm2d/cuda/query_intersect.cu
        farm2d/cuda/query_intersection.cu farm2d/cuda/query_within_polygon.cu)
    target_compile_definitions(tdbase_farm2d PRIVATE USE_GPU)
    target_compile_options(tdbase_farm2d PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda;--expt-relaxed-constexpr>)
    target_link_libraries(tdbase_farm2d PUBLIC CUDA::cuda_driver CUDA::cudart)
endif()
target_link_libraries(tdbase PRIVATE tdbase_farm2d)

function(add_farm2d_target target entrypoint output_name)
    add_executable(${target} ${entrypoint})
    target_link_libraries(${target} PRIVATE tdbase_farm2d)
    set_target_properties(${target} PROPERTIES OUTPUT_NAME ${output_name})
endfunction()
add_farm2d_target(farm2d_intersect_cpu farm2d/query/intersect_cpu.cpp intersect_cpu)
add_farm2d_target(farm2d_intersection_cpu farm2d/query/intersection_cpu.cpp intersection_cpu)
add_farm2d_target(farm2d_within_polygon_cpu farm2d/query/within_polygon_cpu.cpp within_polygon_cpu)
add_farm2d_target(farm2d_dump_polygons farm2d/datagen/dump_polygons.cpp dump_polygons)
add_farm2d_target(farm2d_filter_polygons farm2d/datagen/filter_polygons.cpp filter_polygons)
add_custom_target(farm2d DEPENDS farm2d_intersect_cpu farm2d_intersection_cpu
    farm2d_within_polygon_cpu farm2d_dump_polygons farm2d_filter_polygons)
if(USE_GPU)
    add_farm2d_target(farm2d_intersect_gpu farm2d/query/intersect.cpp intersect)
    add_farm2d_target(farm2d_intersection_gpu farm2d/query/intersection.cpp intersection)
    add_farm2d_target(farm2d_within_polygon_gpu farm2d/query/within_polygon.cpp within_polygon)
    add_dependencies(farm2d farm2d_intersect_gpu farm2d_intersection_gpu farm2d_within_polygon_gpu)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(farm2d_api_test ../tests/farm2d_api.cpp)
    target_link_libraries(farm2d_api_test PRIVATE tdbase_farm2d)
    add_test(NAME farm2d_api COMMAND farm2d_api_test ${CMAKE_CURRENT_BINARY_DIR}/api-validation)
endif()
