add_custom_command(OUTPUT TCATSPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TCATSPhysics.h TCATSPhysicsDict.cxx TCATSPhysics.rootmap libNPCATS.dylib DEPENDS TCATSPhysics.h)
add_custom_command(OUTPUT TCATSDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TCATSData.h TCATSDataDict.cxx TCATSData.rootmap libNPCATS.dylib DEPENDS TCATSData.h)
add_library(NPCATS SHARED TCATSData.cxx TCATSPhysics.cxx TCATSDataDict.cxx TCATSPhysicsDict.cxx TCATSSpectra.cxx)
target_link_libraries(NPCATS ${ROOT_LIBRARIES} NPCore) 
install(FILES TCATSData.h TCATSPhysics.h TCATSSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

