# build DecodeR shared library
add_library(DecodeR SHARED DecodeR.cpp)

# build DecodeD shared library
add_library(DecodeD SHARED DecodeD.cpp)

# build DecodeT shared library
add_library(DecodeT SHARED DecodeT.cpp)

# add the executable
add_executable(online online.cpp)

# link library rules
target_link_libraries(online DecodeR DecodeD DecodeT ${ROOT_LIBRARIES} NPCore NPPhysics NPComptonTelescope)
