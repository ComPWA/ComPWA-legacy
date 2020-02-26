# Check if all requirements are found
if(TARGET EvtGenGenerator)

  add_executable(EvtGenGeneratorTest EvtGenGeneratorTest.cpp)
  target_link_libraries(EvtGenGeneratorTest EvtGenGenerator Core
                        Boost::unit_test_framework)
  set_target_properties(
    EvtGenGeneratorTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                   ${PROJECT_BINARY_DIR}/bin/test/)

  add_test(
    NAME EvtGenGeneratorTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/EvtGenGeneratorTest)
endif()
