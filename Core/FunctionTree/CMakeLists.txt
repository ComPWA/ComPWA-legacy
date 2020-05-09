set(PACKAGE_NAME FunctionTree)

# Define library
set(SOURCE_LIST
    FitParameter.cpp
    Functions.cpp
    FunctionTreeEstimator.cpp
    FunctionTreeIntensity.cpp
    ParameterList.cpp
    TreeNode.cpp
    Value.cpp
)
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
         pstl::ParallelSTL
         Boost::serialization
)

# Install and export
set(HEADER_LIST
    FitParameter.hpp
    Functions.hpp
    FunctionTreeEstimator.hpp
    FunctionTreeIntensity.hpp
    Parameter.hpp
    ParameterList.hpp
    ParObserver.hpp
    TreeNode.hpp
    Value.hpp
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
add_subdirectory(test)
