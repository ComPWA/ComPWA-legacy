set(MODULE_NAME Tools)

# Define library
set(SOURCE_LIST FitFractions.cpp UpdatePTreeParameter.cpp)
add_library(${MODULE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${MODULE_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${MODULE_NAME} PUBLIC Core Integration)

# Install and export
install(
  TARGETS ${MODULE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Define library
set(PACKAGE_NAME Integration)
set(SOURCE_LIST Integration.cpp)
add_library(${PACKAGE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${PACKAGE_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
  ${PACKAGE_NAME}
  PUBLIC Core
         FunctionTree
         Data
         pstl::ParallelSTL
)

# Install and export
set(HEADER_LIST
    FitFractions.hpp
    GoodnessOfFit.hpp
    Integration.hpp
    UpdatePTreeParameter.hpp
)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}
)
install(
  TARGETS ${MODULE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Add unit tests
add_subdirectory(test)

# Add subdirectories
add_subdirectory(Plotting)
add_subdirectory(Adapters)
