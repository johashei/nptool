# build DecodeR shared library
add_library(DecodeR SHARED DecodeR.cpp)

# add the executable
add_executable(online online.cpp)

# link library rules
target_link_libraries(online DecodeR ${ROOT_LIBRARIES} NPCore NPPhysics NPComptonTelescope)
