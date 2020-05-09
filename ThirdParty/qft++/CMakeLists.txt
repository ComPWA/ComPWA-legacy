# Project definition
cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
project(
  qft++
  DESCRIPTION "Numerically calculate quantum field theory expressions."
  # HOMEPAGE_URL "http://www-meg.phys.cmu.edu/williams/qft++/index.php/main_page"
  LANGUAGES CXX
)

# Define library
add_library(${PROJECT_NAME} src/TensorIndex.C src/Tensor.C)
target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
set(QFTPP_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    CACHE INTERNAL ""
)

# Install
set(CMAKE_SHARED_LIBRARY_PREFIX ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/${PROJECT_NAME}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(
  TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
