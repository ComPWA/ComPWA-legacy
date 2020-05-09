set(MODULE_NAME Data)

# Define library
set(SOURCE_LIST DataSet.cpp Generate.cpp)
add_library(${MODULE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${MODULE_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${MODULE_NAME} PUBLIC Core pstl::ParallelSTL)

# Install and export
set(HEADER_LIST DataSet.hpp Generate.hpp)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}
)
install(
  TARGETS ${MODULE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Add subdirectories
add_subdirectory(Ascii)
add_subdirectory(EvtGen)
add_subdirectory(Root)
