add_custom_command(OUTPUT THeliosPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh THeliosPhysics.h THeliosPhysicsDict.cxx THeliosPhysics.rootmap libNPHelios.dylib DEPENDS THeliosPhysics.h)
add_custom_command(OUTPUT THeliosDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh THeliosData.h THeliosDataDict.cxx THeliosData.rootmap libNPHelios.dylib DEPENDS THeliosData.h)
add_library(NPHelios SHARED THeliosSpectra.cxx THeliosData.cxx THeliosPhysics.cxx THeliosDataDict.cxx THeliosPhysicsDict.cxx )
target_link_libraries(NPHelios ${ROOT_LIBRARIES} NPCore) 
install(FILES THeliosData.h THeliosPhysics.h THeliosSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

