add_custom_command(OUTPUT THiraPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh THiraPhysics.h THiraPhysicsDict.cxx THiraPhysics.rootmap libNPHira.dylib DEPENDS THiraPhysics.h)
add_custom_command(OUTPUT THiraDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh THiraData.h THiraDataDict.cxx THiraData.rootmap libNPHira.dylib DEPENDS THiraData.h)
add_library(NPHira SHARED THiraData.cxx THiraPhysics.cxx THiraDataDict.cxx THiraPhysicsDict.cxx )
target_link_libraries(NPHira ${ROOT_LIBRARIES} NPCore) 
install(FILES THiraData.h THiraPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

