add_library(NPSIss SHARED  IssArray.cc)
target_link_libraries(NPSIss NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPIss)
