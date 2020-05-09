set(SOURCE_LIST GFMinIndividual.cpp GenevaIF.cpp)
set(HEADER_LIST
    GFMinIndividual.hpp
    GenevaIF.hpp
    GenevaResult.hpp
)

add_library(GenevaIF ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(
  GenevaIF
  PUBLIC Core Boost::serialization
  PRIVATE ${GENEVA_LIBRARIES}
          Boost::program_options
          Boost::filesystem
)

target_include_directories(GenevaIF PUBLIC ${GENEVA_INCLUDE_DIR})

install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Optimizer/Geneva
)
install(
  TARGETS GenevaIF
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY config DESTINATION share/Geneva)

add_subdirectory(test)
