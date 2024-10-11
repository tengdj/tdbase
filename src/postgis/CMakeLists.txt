cmake_minimum_required(VERSION 3.16)
project(POSTGIS_TEST)

# find_package(pqxx REQUIRED)
find_package(PostgreSQL REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
# set(Boost_INCLUDE_DIRS "/usr/include/boost")
# set(Boost_LIBRARIES )
set(PQXX_INCLUDE_DIRS "/usr/include/pqxx")
set(PQXX_LIBRARIES "/usr/lib/x86_64-linux-gnu/libpqxx.so")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -Wno-unused-result")
SET (BOOST_ROOT "${BOOST_DIR}")
message("boost lib: ${Boost_LIBRARIES}")
message("boost inc: ${Boost_INCLUDE_DIR}")

# add_executable(postgis_test src/main.cpp src/config.h)
# target_include_directories(postgis_test PRIVATE ${PostgreSQL_INCLUDE_DIRS} ${PQXX_INCLUDE_DIRS})
# target_link_libraries(postgis_test PRIVATE pqxx ${PostgreSQL_LIBRARIES} ${PQXX_LIBRARIES})
# target_link_libraries(postgis_test pqxx)

# add_executable(within src/within.cpp src/config.h)
# target_include_directories(within PRIVATE ${PostgreSQL_INCLUDE_DIRS} ${PQXX_INCLUDE_DIRS})
# target_link_libraries(within PRIVATE ${PostgreSQL_LIBRARIES} ${PQXX_LIBRARIES})

# add_executable(within_origin src/within_origin.cpp src/config.h)
# target_include_directories(within_origin PRIVATE ${PostgreSQL_INCLUDE_DIRS} ${PQXX_INCLUDE_DIRS})
# target_link_libraries(within_origin PRIVATE ${PostgreSQL_LIBRARIES} ${PQXX_LIBRARIES})

# add_executable(nn src/nn.cpp src/config.h)
# target_include_directories(nn PRIVATE ${PostgreSQL_INCLUDE_DIRS} ${PQXX_INCLUDE_DIRS})
# target_link_libraries(nn PRIVATE ${PostgreSQL_LIBRARIES} ${PQXX_LIBRARIES})

# add_executable(nn_origin src/nn_origin.cpp src/config.h)
# target_include_directories(nn_origin PRIVATE ${PostgreSQL_INCLUDE_DIRS} ${PQXX_INCLUDE_DIRS})
# target_link_libraries(nn_origin PRIVATE ${PostgreSQL_LIBRARIES} ${PQXX_LIBRARIES})

add_executable(main src/main.cpp src/config.h)
target_include_directories(main PUBLIC ${PostgreSQL_INCLUDE_DIRS} ${PQXX_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
target_link_libraries(main PUBLIC ${PostgreSQL_LIBRARIES}  ${PQXX_LIBRARIES} ${Boost_LIBRARIES})