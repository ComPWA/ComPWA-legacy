# Check if all requirements are found
if(TARGET Tools)

  add_executable(UpdatePTreeParameterTest UpdatePTreeParameterTest.cpp)
  target_link_libraries(UpdatePTreeParameterTest Tools
                        Boost::unit_test_framework)
  set_target_properties(
    UpdatePTreeParameterTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                        ${PROJECT_BINARY_DIR}/bin/test/)

  add_test(
    NAME UpdatePTreeParameterTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/UpdatePTreeParameterTest)

  add_executable(IntegrationTest IntegrationTest.cpp)

  target_link_libraries(IntegrationTest Data RootData Tools
                        Boost::unit_test_framework)
  set_target_properties(
    IntegrationTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                               ${PROJECT_BINARY_DIR}/bin/test/)
  add_test(
    NAME IntegrationTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/IntegrationTest)

  add_executable(GoodnessOfFitTest GoodnessOfFitTest.cpp)
  target_link_libraries(GoodnessOfFitTest Tools Boost::unit_test_framework)
  set_target_properties(
    GoodnessOfFitTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                 ${PROJECT_BINARY_DIR}/bin/test/)
  add_test(
    NAME GoodnessOfFitTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/GoodnessOfFitTest)
endif()
