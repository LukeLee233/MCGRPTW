cmake_minimum_required(VERSION 3.15)
project(MCGRP_debug)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(MCGRP_debug LANGUAGES CXX)

string(TIMESTAMP BUILDDATE "%Y-%m-%d %H:%M:%S")
message(STATUS "make date: ${BUILDDATE}")
message(STATUS "start running cmake... ")
set(CMAKE_CXX_FLAGS -pthread)

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
find_package(GSL)

message(STATUS "gsl_INCLUDE_DIRS: ${GSL_INCLUDE_DIR}")
message(STATUS "gsl_LIBRARIES: ${GSL_LIBRARIES}")
message(STATUS "gsl_VERSION: ${GSL_VERSION}")

message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
message(STATUS "Boost_VERSION: ${Boost_VERSION}")

#add_definitions(-D DEBUG)

# disable compile warning flag
add_compile_options(-w)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIR})

aux_source_directory(. MCGRP_src)

list(REMOVE_ITEM MCGRP_src ./Memetic.h ./Memetic.cpp)

add_executable(MCGRP_debug ${MCGRP_src})
target_link_libraries(MCGRP_debug gsl gslcblas ${Boost_LIBRARIES})
