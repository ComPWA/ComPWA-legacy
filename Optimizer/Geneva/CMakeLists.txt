set(lib_srcs GFMinIndividual.cpp GenevaIF.cpp)
set(lib_headers GFMinIndividual.hpp GenevaIF.hpp GenevaResult.hpp)

add_library(GenevaIF ${lib_srcs} ${lib_headers})

target_link_libraries(
  GenevaIF
  PUBLIC Core Boost::serialization
  PRIVATE ${GENEVA_LIBRARIES} Boost::program_options Boost::filesystem)

target_include_directories(GenevaIF PUBLIC ${GENEVA_INCLUDE_DIR})

install(FILES ${lib_headers} DESTINATION include/Optimizer/Geneva)
install(
  TARGETS GenevaIF
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
install(DIRECTORY config DESTINATION share/Geneva)

add_subdirectory(test)
