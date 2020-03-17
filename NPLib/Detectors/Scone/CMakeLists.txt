add_custom_command(OUTPUT TSconePhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TSconePhysics.h TSconePhysicsDict.cxx TSconePhysics.rootmap libNPScone.dylib DEPENDS TSconePhysics.h)
add_custom_command(OUTPUT TSconeDataDict.cxx COMMAND ../../scripts/build_dict.sh TSconeData.h TSconeDataDict.cxx TSconeData.rootmap libNPScone.dylib DEPENDS TSconeData.h)
add_library(NPScone SHARED TSconeSpectra.cxx TSconeData.cxx TSconePhysics.cxx TSconeDataDict.cxx TSconePhysicsDict.cxx )
target_link_libraries(NPScone ${ROOT_LIBRARIES} NPCore) 
install(FILES TSconeData.h TSconePhysics.h TSconeSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

