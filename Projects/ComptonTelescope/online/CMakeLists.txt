cmake_minimum_required (VERSION 2.6)
project (CCAM)
set (CCAM_VERSION_MAJOR 0)
set (CCAM_VERSION_MINOR 1)
set(CMAKE_BUILD_TYPE Release)

# support for c++11
#set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11 -m64")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -m64")

# set directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_BINARY_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/bin )

# add include directory
include_directories(src)

# add source directory
add_subdirectory(src)

# Link executable against the local lib
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

