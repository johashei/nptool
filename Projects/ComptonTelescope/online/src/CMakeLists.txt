# build DecodeR shared library
add_library(DecodeR SHARED DecodeR.cpp)

# build DecodeR shared library
add_library(DecodeD SHARED DecodeD.cpp)

# add the executable
add_executable(online online.cpp)

# link library rules
target_link_libraries(online DecodeR DecodeD ${ROOT_LIBRARIES} NPCore NPPhysics NPComptonTelescope)
