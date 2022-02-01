add_library(NPSSofTofW SHARED  SofTofW.cc)
add_library(NPSSofTwim SHARED  SofTwim.cc)
target_link_libraries(NPSSofTofW NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPSofia)
target_link_libraries(NPSSofTwim NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPSofia)
