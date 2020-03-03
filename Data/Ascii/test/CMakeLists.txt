if(TARGET Ascii)

  function(copy_test_file LOCAL_FILE)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_FILE}
                   ${PROJECT_BINARY_DIR}/bin/test/Data_${LOCAL_FILE} COPYONLY)
  endfunction()

  # AsciiDataIOTest
  add_executable(Data_AsciiDataIOTest AsciiDataIOTest.cpp)
  target_link_libraries(Data_AsciiDataIOTest Ascii EvtGenGenerator
                        Boost::unit_test_framework)
  set_target_properties(
    Data_AsciiDataIOTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                    ${PROJECT_BINARY_DIR}/bin/test/)
  add_test(
    NAME Data_AsciiDataIOTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiDataIOTest)
  copy_test_file(AsciiDataIOTest-CorrectWeightedEmom.dat)
  copy_test_file(AsciiDataIOTest-CorrectWeightedMomE.dat)
  copy_test_file(AsciiDataIOTest-CorrectUnweighted.dat)

  # AsciiHeaderIOTest
  add_executable(Data_AsciiHeaderIOTest AsciiHeaderIOTest.cpp)
  target_link_libraries(Data_AsciiHeaderIOTest Ascii EvtGenGenerator
                        Boost::unit_test_framework)
  set_target_properties(
    Data_AsciiHeaderIOTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                      ${PROJECT_BINARY_DIR}/bin/test/)
  add_test(
    NAME Data_AsciiHeaderIOTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiHeaderIOTest)
  copy_test_file(AsciiHeaderIOTest-YAMLcorrect.dat)

endif()
