add_custom_command(OUTPUT TBigRIPSPPACDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPPACData.h TBigRIPSPPACDataDict.cxx TBigRIPSPPACData.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPPACData.h)

add_custom_command(OUTPUT TBigRIPSPPACPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPPACPhysics.h TBigRIPSPPACPhysicsDict.cxx TBigRIPSPPACPhysics.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPPACPhysics.h)

add_library(NPBigRIPS SHARED TBigRIPSPPACData.cxx TBigRIPSPPACDataDict.cxx TBigRIPSPPACPhysics.cxx TBigRIPSPPACPhysicsDict.cxx)

target_link_libraries(NPBigRIPS ${ROOT_LIBRARIES} NPCore NPTrackReconstruction) 
install(FILES TBigRIPSPPACData.h TBigRIPSPPACPhysics.h  DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

