set(PACKAGE_NAME EvtGenIF)

# Define library
set(SOURCE_LIST DalitzKinematics.cpp EvtGenIF.cpp)
add_library(${PACKAGE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${PACKAGE_NAME} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
  ${PACKAGE_NAME}
  PUBLIC Core
         FunctionTree
         Physics
         EvtGen
         qft++
)

# Install and export
set(HEADER_LIST DalitzKinematics.hpp EvtGenIF.hpp)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/EvtGen
)
install(
  TARGETS ${PACKAGE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
