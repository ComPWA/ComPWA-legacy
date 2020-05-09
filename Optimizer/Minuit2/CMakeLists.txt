set(PACKAGE_NAME Minuit2IF)

# Define library
set(SOURCE_LIST MinuitIF.cpp MinuitResult.cpp)
add_library(${PACKAGE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${PACKAGE_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${ROOT_INCLUDE_DIR}
)
target_link_libraries(
  ${PACKAGE_NAME}
  PUBLIC Core
  PRIVATE Minuit2::Minuit2
)

# Install and export
set(HEADER_LIST
    MinuitFcn.hpp
    MinuitIF.hpp
    MinuitResult.hpp
)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/Minuit2
)
install(
  TARGETS ${PACKAGE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Add unit tests
add_subdirectory(test)
