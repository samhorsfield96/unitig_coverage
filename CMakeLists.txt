cmake_minimum_required(VERSION 3.10)
project(Unitig_coverage)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -march=native")

set(CMAKE_VERBOSE_MAKEFILE ON)

# find pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# check for conda environment
IF( DEFINED ENV{CONDA_PREFIX} )
    # set conda directory
    SET(CONDA_DIR "$ENV{CONDA_PREFIX}")
ENDIF()

#Find include dirs for z and bifrost
FIND_PATH(z_INCLUDE_DIR NAMES zlib.h HINTS ${CONDA_DIR}/include)
FIND_PATH(bifrost_INCLUDE_DIR bifrost/CompactedDBG.hpp HINTS ${CONDA_DIR}/include)

# Find libraries for z and bifrost
FIND_LIBRARY(z_LIBRARY NAMES z HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_LIBRARY(bifrost_LIBRARY NAMES bifrost HINTS ${CONDA_DIR}/lib REQUIRED)

# Set package found as true
FIND_PACKAGE_HANDLE_STANDARD_ARGS(z DEFAULT_MSG z_INCLUDE_DIR z_LIBRARY)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(bifrost DEFAULT_MSG bifrost_INCLUDE_DIR bifrost_LIBRARY)

# Set dirs if multiple include/libraries found
SET( z_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( z_LIBRARIES ${z_LIBRARY} )
SET( bifrost_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( bifrost_LIBRARIES ${bifrost_LIBRARY} )

add_executable(Unitig_coverage main.cpp)

# link libraries
target_link_libraries(Unitig_coverage PRIVATE ${bifrost_LIBRARIES} ${z_LIBRARIES} Threads::Threads)



