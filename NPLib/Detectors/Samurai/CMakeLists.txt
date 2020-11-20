add_custom_command(OUTPUT TSamuraiFDC2DataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC2Data.h TSamuraiFDC2DataDict.cxx TSamuraiFDC2Data.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC2Data.h)

add_custom_command(OUTPUT TSamuraiFDC2PhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC2Physics.h TSamuraiFDC2PhysicsDict.cxx TSamuraiFDC2Physics.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC2Physics.h)

add_custom_command(OUTPUT TSamuraiFDC0DataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC0Data.h TSamuraiFDC0DataDict.cxx TSamuraiFDC0Data.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC0Data.h)

add_custom_command(OUTPUT TSamuraiFDC0PhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC0Physics.h TSamuraiFDC0PhysicsDict.cxx TSamuraiFDC0Physics.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC0Physics.h)


add_library(NPSamurai SHARED TSamuraiFDC2Data.cxx TSamuraiFDC2DataDict.cxx TSamuraiFDC2Physics.cxx TSamuraiFDC2PhysicsDict.cxx TSamuraiFDC0Data.cxx TSamuraiFDC0DataDict.cxx TSamuraiFDC0Physics.cxx TSamuraiFDC0PhysicsDict.cxx)

target_link_libraries(NPSamurai ${ROOT_LIBRARIES} NPCore NPTrackReconstruction) 
install(FILES TSamuraiFDC2Data.h TSamuraiFDC2Physics.h  TSamuraiFDC0Data.h TSamuraiFDC0Physics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

