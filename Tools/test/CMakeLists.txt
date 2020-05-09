# UpdatePTreeParameterTest
set(TEST_NAME UpdatePTreeParameterTest)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${TEST_NAME} PRIVATE Tools Boost::unit_test_framework)
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
)
add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
)

# IntegrationTest
if(TARGET RootData)
  set(TEST_NAME IntegrationTest)
  add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
  target_include_directories(
    ${TEST_NAME} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_link_libraries(
    ${TEST_NAME}
    PRIVATE Data
            RootData
            Tools
            Boost::unit_test_framework
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
endif()

# GoodnessOfFitTest
set(TEST_NAME GoodnessOfFitTest)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(${TEST_NAME} PRIVATE Tools Boost::unit_test_framework)
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
)
add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
)
