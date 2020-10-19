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
link_directories(${NPLIB}/lib)

# build DecodeR shared library
add_library(DecodeR SHARED DecodeR.cpp)

# add the executable
add_executable(online online.cpp)

# link library rules
#target_link_libraries(online DecodeR ${ROOT_LIBRARIES} NPCore NPPhysics NPComptonTelescope)
target_link_libraries(online ${ROOT_LIBRARIES} NPCore NPPhysics)
