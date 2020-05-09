# ModelFunctionTreeTests
set(TEST_NAME ModelFunctionTreeTests)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${ROOT_INCLUDE_DIR}
         ${Boost_INCLUDE_DIR}
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Boost::unit_test_framework
          Core
          Physics
          RootData
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

# HelicityAngle
set(TEST_NAME HelicityAngle)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${Boost_INCLUDE_DIR}
         ${QFTPP_INCLUDE_DIR}
         ${ROOT_INCLUDE_DIR}
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Boost::unit_test_framework
          Core
          Physics
          RootData
          qft++
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

# HelicityAnglesCorrectnessTest
set(TEST_NAME HelicityAnglesCorrectnessTest)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${Boost_INCLUDE_DIR}
         ${QFTPP_INCLUDE_DIR}
         ${ROOT_INCLUDE_DIR}
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Boost::unit_test_framework
          Core
          Physics
          RootData
          qft++
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

# HelicityKinematicsTests
set(TEST_NAME HelicityKinematicsTests)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${Boost_INCLUDE_DIR}
         ${ROOT_INCLUDE_DIR}
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Boost::unit_test_framework
          Core
          Physics
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
