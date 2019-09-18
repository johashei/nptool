include(CheckLanguage)

add_custom_command(OUTPUT TMinosPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TMinosPhysics.h TMinosPhysicsDict.cxx TMinosPhysics.rootmap libNPMinos.dylib DEPENDS TMinosPhysics.h)
add_custom_command(OUTPUT TMinosDataDict.cxx COMMAND ../../scripts/build_dict.sh TMinosData.h TMinosDataDict.cxx TMinosData.rootmap libNPMinos.dylib DEPENDS TMinosData.h)
add_custom_command(OUTPUT TMinosClustDict.cxx COMMAND ../../scripts/build_dict.sh TMinosClust.h TMinosClustDict.cxx TMinosClust.rootmap libNPMinos.dylib DEPENDS TMinosClust.h)
add_custom_command(OUTPUT TMinosResultDict.cxx COMMAND ../../scripts/build_dict.sh TMinosResult.h TMinosResultDict.cxx TMinosResult.rootmap libNPMinos.dylib DEPENDS TMinosResult.h)

## Check if MINUIT2 is installed along with ROOT
find_library(libMinuit_FOUND NAMES Minuit HINTS "$ENV{ROOTSYS}/lib")
if(libMinuit_FOUND)
  message(STATUS "Minuit support enabled for Minos.")
  add_definitions(-DHAVE_MINUIT)
else()
  message(STATUS "Minuit support disabled for Minos. Not found in $ENV{ROOTSYS}/lib")
endif()

add_library(NPMinos SHARED TMinosClust.cxx TMinosData.cxx TMinosPhysics.cxx TMinosResult TMinosDataDict.cxx TMinosPhysicsDict.cxx TMinosClustDict.cxx TMinosResultDict.cxx) 

if(libMinuit_FOUND)
## Link to Minuit2 library
  target_link_libraries(NPMinos ${ROOT_LIBRARIES} NPCore NPTrackReconstruction Minuit)
else()
## No Minuit2 library, skip linking     
  target_link_libraries(NPMinos ${ROOT_LIBRARIES} NPCore  NPTrackReconstruction)
endif()

install(FILES TMinosClust.h TMinosData.h TMinosPhysics.h TMinosResult.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

