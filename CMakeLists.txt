cmake_minimum_required(VERSION 3.10)
project(Unitig_coverage)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -march=native")

set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories("src/include")

find_library(pthread REQUIRED)
find_library(z REQUIRED)
find_library(bifrost REQUIRED)

add_executable(Unitig_coverage main.cpp)

target_link_libraries(Unitig_coverage PRIVATE pthread z bifrost)
