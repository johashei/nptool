add_custom_command(OUTPUT TSAMURAIPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TSAMURAIPhysics.h TSAMURAIPhysicsDict.cxx TSAMURAIPhysics.rootmap libNPSAMURAI.dylib DEPENDS TSAMURAIPhysics.h)
add_custom_command(OUTPUT TSAMURAIDataDict.cxx COMMAND ../../scripts/build_dict.sh TSAMURAIData.h TSAMURAIDataDict.cxx TSAMURAIData.rootmap libNPSAMURAI.dylib DEPENDS TSAMURAIData.h)
add_library(NPSAMURAI SHARED TSAMURAISpectra.cxx TSAMURAIData.cxx TSAMURAIPhysics.cxx TSAMURAIDataDict.cxx TSAMURAIPhysicsDict.cxx )
target_link_libraries(NPSAMURAI ${ROOT_LIBRARIES} NPCore) 
install(FILES TSAMURAIData.h TSAMURAIPhysics.h TSAMURAISpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

