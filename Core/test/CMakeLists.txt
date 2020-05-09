set(TEST_NAME "CoreTest_FourMomentumTest")

add_executable(${TEST_NAME} FourMomentumTest.cpp)

target_link_libraries(
  ${TEST_NAME}
  Core
  Threads::Threads
  Boost::unit_test_framework
)

target_include_directories(${TEST_NAME} PRIVATE Boost::unit_test_framework)

# Move testing binaries into a testBin directory
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
  COMMAND ${PROJECT_BINARY_DIR}/bin/test/${TEST_NAME}
)

set(TEST_NAME "CoreTest_PropertiesTest")

add_executable(${TEST_NAME} PropertiesTest.cpp)

target_link_libraries(
  ${TEST_NAME}
  Core
  Threads::Threads
  Boost::unit_test_framework
)

target_include_directories(${TEST_NAME} PRIVATE Boost::unit_test_framework)

# Move testing binaries into a testBin directory
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
  COMMAND ${PROJECT_BINARY_DIR}/bin/test/${TEST_NAME}
)
