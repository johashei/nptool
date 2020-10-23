cmake_minimum_required (VERSION 2.6)
project (CCAM)
set (CCAM_VERSION_MAJOR 0)
set (CCAM_VERSION_MINOR 1)
set(CMAKE_BUILD_TYPE Release)

# ROOT stuff                                                                    
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})                                    
find_package(ROOT REQUIRED)                                                     
include_directories(${ROOT_INCLUDE_DIRS})                                       
                                                                                
# nptool stuff                                                                  
find_package(NPLib)                                                             
include(${NPLib_USE_FILE})                                                      
message("  found: ${NPLIB}")                                                    
message("  include directory: ${NPLib_INCLUDE_DIRS}")                           
message("  detector list (empty if all): ${NPLib_DETECTOR_LIST}")               
message("  library directory: ${NPLib_LIB_DIRS}")                               
link_directories(${NPLIB}/lib)                                                  

# support for c++11
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# set directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_BINARY_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/bin )

# add include directory
include_directories(src)

# add source directory
add_subdirectory(src)

# Link executable against the local lib
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

