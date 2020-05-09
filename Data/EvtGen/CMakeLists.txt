set(PACKAGE_NAME EvtGenGenerator)

# Define library
set(SOURCE_LIST EvtGenGenerator.cpp)
add_library(${PACKAGE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${PACKAGE_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${Boost_INCLUDE_DIR}
)
target_link_libraries(
  ${PACKAGE_NAME}
  PUBLIC Core
         Data
         EvtGen
         Physics
)

# Install and export
set(HEADER_LIST EvtGenGenerator.hpp)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/EvtGen
)
install(
  TARGETS ${PACKAGE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Add unit tests
add_subdirectory(test)
