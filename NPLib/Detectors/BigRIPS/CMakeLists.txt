add_custom_command(OUTPUT TBigRIPSPPACDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPPACData.h TBigRIPSPPACDataDict.cxx TBigRIPSPPACData.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPPACData.h)

add_custom_command(OUTPUT TBigRIPSPPACPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPPACPhysics.h TBigRIPSPPACPhysicsDict.cxx TBigRIPSPPACPhysics.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPPACPhysics.h)

add_custom_command(OUTPUT TBigRIPSPlasticDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPlasticData.h TBigRIPSPlasticDataDict.cxx TBigRIPSPlasticData.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPlasticData.h)

add_custom_command(OUTPUT TBigRIPSPlasticPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPlasticPhysics.h TBigRIPSPlasticPhysicsDict.cxx TBigRIPSPlasticPhysics.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPlasticPhysics.h)

add_library(NPBigRIPS SHARED TBigRIPSPPACData.cxx TBigRIPSPPACDataDict.cxx TBigRIPSPPACPhysics.cxx TBigRIPSPPACPhysicsDict.cxx TBigRIPSPlasticData.cxx TBigRIPSPlasticDataDict.cxx TBigRIPSPlasticPhysics.cxx TBigRIPSPlasticPhysicsDict.cxx )

target_link_libraries(NPBigRIPS ${ROOT_LIBRARIES} NPCore NPTrackReconstruction) 
install(FILES TBigRIPSPPACData.h TBigRIPSPPACPhysics.h TBigRIPSPlasticData.h TBigRIPSPPACPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

