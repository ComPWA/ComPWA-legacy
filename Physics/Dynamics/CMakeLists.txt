set(PACKAGE_NAME Dynamics)

# Define library
set(SOURCE_LIST
    Flatte.cpp
    FormFactor.cpp
    RelativisticBreitWigner.cpp
    Utils/Faddeeva.cc
    Voigtian.cpp
)
add_library(${PACKAGE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${PACKAGE_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${PACKAGE_NAME} PUBLIC Core FunctionTree)

# Install and export
set(HEADER_LIST
    Coupling.hpp
    Flatte.hpp
    FormFactor.hpp
    RelativisticBreitWigner.hpp
    Voigtian.hpp
)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/${PACKAGE_NAME}
)
install(
  TARGETS ${PACKAGE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Add unit tests
# add_subdirectory(test)
