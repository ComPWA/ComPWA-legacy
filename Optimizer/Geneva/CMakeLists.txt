set(SOURCE_LIST GFMinIndividual.cpp GenevaIF.cpp)
set(HEADER_LIST
    GFMinIndividual.hpp
    GenevaIF.hpp
    GenevaResult.hpp
)

add_library(GenevaIF ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(
  GenevaIF
  PUBLIC Boost::serialization Core
  PRIVATE ${GENEVA_LIBRARIES}
          Boost::filesystem
          Boost::program_options
)

target_include_directories(GenevaIF PUBLIC ${GENEVA_INCLUDE_DIR})

install(FILES ${HEADER_LIST} DESTINATION include/Optimizer/Geneva)
install(
  TARGETS GenevaIF
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY config DESTINATION share/Geneva)

add_subdirectory(test)
