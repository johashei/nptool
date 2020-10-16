# ROOT stuff
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#find_package(ROOT REQUIRED)
#include_directories(${ROOT_INCLUDE_DIRS})

# build libStarClusterMC shared library
add_library(DecodeR SHARED DecodeR.cpp)

# add the executable
add_executable(online online.cpp)

# link library rules
target_link_libraries(online DecodeR)
