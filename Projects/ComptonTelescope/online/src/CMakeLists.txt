# ROOT stuff
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

# nptool stuff
#list(APPEND CMAKE_PREFIX_PATH $ENV{NPTOOL})
find_package(NPLib)
include(${NPLib_USE_FILE})
message("  found: ${NPLIB}")
message("  include directory: ${NPLib_INCLUDE_DIRS}")
message("  detector list (empty if all): ${NPLib_DETECTOR_LIST}")
message("  library directory: ${NPLib_LIB_DIRS}")

#set(NPTOOL "$ENV{NPTOOL}")                                                      
#set(NPLIB "${NPTOOL}/NPLib")                                                    
#set(NPTOOL_INCLUDE_DIR "${NPLIB}/include")                                      
set(NPTOOL_LIB_DIR "${NPLIB}/lib")                                              
#include_directories("${NPTOOL_INCLUDE_DIR}")
#set(NPLib_LIBRARIES "-L${NPTOOL_LIB_DIR} -lNPCore -lNPPhysics -lNPComptonTelescope")
set(NPLib_LIBRARIES "-L${NPTOOL_LIB_DIR} -lNPCore -lNPPhysics")

# build libStarClusterMC shared library
add_library(DecodeR SHARED DecodeR.cpp)

# add the executable
add_executable(online online.cpp)

# link library rules
#target_link_libraries(online DecodeR ${ROOT_LIBRARIES} ${NPLib_LIBRARIES})
target_link_libraries(online ${ROOT_LIBRARIES} ${NPLib_LIBRARIES})
#target_link_libraries(online DecodeR ${ROOT_LIBRARIES} NPCore NPPhysics NPComptonTelescope)
