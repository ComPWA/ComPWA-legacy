# MinuitTestApp
set(TEST_NAME MinuitTestApp)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Boost::unit_test_framework
          Core
          Minuit2IF
)
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
)
add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
)

# FitResultSerializationTest
set(TEST_NAME FitResultSerializationTest)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Boost::unit_test_framework
          Core
          Minuit2IF
)
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
)
add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
)
