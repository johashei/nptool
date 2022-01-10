add_library(NPSSofTofW SHARED  SofTofW.cc)
target_link_libraries(NPSSofTofW NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPSofia)
