add_library(NPSSamurai SHARED  Samurai.cc SamuraiFieldPropagation.cc)
target_link_libraries(NPSSamurai NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPSamurai)
