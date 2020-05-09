# Helper function
function(copy_test_file LOCAL_FILE)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_FILE}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/Data_${LOCAL_FILE}
    COPYONLY
  )
endfunction()

# AsciiDataIOTest
set(TEST_NAME AsciiDataIOTest)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${Boost_INCLUDE_DIR}
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Ascii
          Boost::unit_test_framework
          EvtGenGenerator
)
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
)
add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test
  COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
)
copy_test_file(AsciiDataIOTest-CorrectWeightedEmom.dat)
copy_test_file(AsciiDataIOTest-CorrectWeightedMomE.dat)
copy_test_file(AsciiDataIOTest-CorrectUnweighted.dat)

# AsciiHeaderIOTest
set(TEST_NAME AsciiHeaderIOTest)
add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
target_include_directories(
  ${TEST_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${Boost_INCLUDE_DIR}
)
target_link_libraries(
  ${TEST_NAME}
  PRIVATE Ascii
          Boost::unit_test_framework
          EvtGenGenerator
)
set_target_properties(
  ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
)
add_test(
  NAME ${TEST_NAME}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test
  COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
)
copy_test_file(AsciiHeaderIOTest-YAMLcorrect.dat)
