add_custom_command(OUTPUT TChateauCristalPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TChateauCristalPhysics.h TChateauCristalPhysicsDict.cxx TChateauCristalPhysics.rootmap libNPChateauCristal.dylib DEPENDS TChateauCristalPhysics.h)
add_custom_command(OUTPUT TChateauCristalDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TChateauCristalData.h TChateauCristalDataDict.cxx TChateauCristalData.rootmap libNPChateauCristal.dylib DEPENDS TChateauCristalData.h)
add_library(NPChateauCristal SHARED TChateauCristalData.cxx TChateauCristalPhysics.cxx TChateauCristalDataDict.cxx TChateauCristalPhysicsDict.cxx )
target_link_libraries(NPChateauCristal ${ROOT_LIBRARIES} NPCore) 
install(FILES TChateauCristalData.h TChateauCristalPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

