add_custom_command(OUTPUT TParisPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TParisPhysics.h TParisPhysicsDict.cxx TParisPhysics.rootmap libNPParis.dylib DEPENDS TParisPhysics.h)
add_custom_command(OUTPUT TParisDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TParisData.h TParisDataDict.cxx TParisData.rootmap libNPParis.dylib DEPENDS TParisData.h)
add_library(NPParis SHARED TParisData.cxx TParisPhysics.cxx TParisDataDict.cxx TParisPhysicsDict.cxx )
target_link_libraries(NPParis ${ROOT_LIBRARIES} NPCore) 
install(FILES TParisData.h TParisPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

