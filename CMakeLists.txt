cmake_minimum_required(VERSION 3.15)
project(MCGRPTW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(MCGRPTW LANGUAGES CXX)



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

add_definitions(-D DEBUG)

# disable compile warning flag
add_compile_options(-w)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIR})


aux_source_directory(. src)

list(REMOVE_ITEM src ./Memetic.h ./Memetic.cpp ./similarity.h ./similarity.cpp)

add_executable(MCGRPTW ${src})
target_link_libraries(MCGRPTW
        gsl
        gslcblas
        ${Boost_LIBRARIES})
