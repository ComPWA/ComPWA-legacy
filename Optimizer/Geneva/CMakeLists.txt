set(PACKAGE_NAME GenevaIF)

# Define library
set(SOURCE_LIST GFMinIndividual.cpp GenevaIF.cpp)
add_library(${PACKAGE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${PACKAGE_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${Boost_INCLUDE_DIR} ${GENEVA_INCLUDE_DIR}
)
target_link_libraries(
  ${PACKAGE_NAME}
  PUBLIC Core Boost::serialization
  PRIVATE ${GENEVA_LIBRARIES}
          Boost::program_options
          Boost::filesystem
)

# Install and export
set(HEADER_LIST
    GFMinIndividual.hpp
    GenevaIF.hpp
    GenevaResult.hpp
)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/Geneva
)
install(
  TARGETS ${PACKAGE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY config DESTINATION share/Geneva)

# Add unit tests
add_subdirectory(test)
