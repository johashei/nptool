add_custom_command(OUTPUT TSofSciDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofSciData.h TSofSciDataDict.cxx TSofSciData.rootmap libNPSofia.dylib DEPENDS TSofSciData.h)

add_custom_command(OUTPUT TSofTofWDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTofWData.h TSofTofWDataDict.cxx TSofTofWData.rootmap libNPSofia.dylib DEPENDS TSofTofWData.h)
add_custom_command(OUTPUT TSofTofWPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTofWPhysics.h TSofTofWPhysicsDict.cxx TSofTofWPhysics.rootmap libNPSofia.dylib DEPENDS TSofTofWPhysics.h)

add_custom_command(OUTPUT TSofMwpcDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofMwpcData.h TSofMwpcDataDict.cxx TSofMwpcData.rootmap libNPSofia.dylib DEPENDS TSofMwpcData.h)

add_custom_command(OUTPUT TSofAtDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofAtData.h TSofAtDataDict.cxx TSofAtData.rootmap libNPSofia.dylib DEPENDS TSofAtData.h)

add_custom_command(OUTPUT TSofTrimDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTrimData.h TSofTrimDataDict.cxx TSofTrimData.rootmap libNPSofia.dylib DEPENDS TSofTrimData.h)
add_custom_command(OUTPUT TSofTrimPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTrimPhysics.h TSofTrimPhysicsDict.cxx TSofTrimPhysics.rootmap libNPSofia.dylib DEPENDS TSofTrimPhysics.h)

add_custom_command(OUTPUT TSofTwimDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTwimData.h TSofTwimDataDict.cxx TSofTwimData.rootmap libNPSofia.dylib DEPENDS TSofTwimData.h)

add_library(NPSofia SHARED TSofSciData.cxx TSofSciDataDict.cxx TSofMwpcData.cxx TSofMwpcDataDict.cxx TSofAtData.cxx TSofAtDataDict.cxx TSofTrimData.cxx TSofTrimDataDict.cxx TSofTrimPhysics.cxx TSofTrimPhysicsDict.cxx TSofTwimData.cxx TSofTwimDataDict.cxx TSofTofWData.cxx TSofTofWDataDict.cxx TSofTofWPhysics.cxx TSofTofWPhysicsDict.cxx)

target_link_libraries(NPSofia ${ROOT_LIBRARIES} NPCore NPPhysics)

install(FILES TSofSciData.h TSofMwpcData.h TSofAtData.h TSofTrimData.h TSofTrimPhysics.h TSofTwimData.h TSofTofWData.h TSofTofWPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

