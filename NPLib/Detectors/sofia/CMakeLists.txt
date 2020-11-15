add_custom_command(OUTPUT TSofiaPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofiaPhysics.h TSofiaPhysicsDict.cxx TSofiaPhysics.rootmap libNPSofia.dylib DEPENDS TSofiaPhysics.h)
add_custom_command(OUTPUT TSofiaDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofiaData.h TSofiaDataDict.cxx TSofiaData.rootmap libNPSofia.dylib DEPENDS TSofiaData.h)
add_library(NPSofia SHARED TSofiaSpectra.cxx TSofiaData.cxx TSofiaPhysics.cxx TSofiaDataDict.cxx TSofiaPhysicsDict.cxx )
target_link_libraries(NPSofia ${ROOT_LIBRARIES} NPCore) 
install(FILES TSofiaData.h TSofiaPhysics.h TSofiaSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

