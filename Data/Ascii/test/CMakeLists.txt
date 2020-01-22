if(TARGET Ascii)

  add_executable(Data_AsciiDataIOTest AsciiDataIOTest.cpp)

  target_link_libraries(Data_AsciiDataIOTest
    Ascii EvtGenGenerator
    Boost::unit_test_framework
  )
  set_target_properties(Data_AsciiDataIOTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(NAME Data_AsciiDataIOTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiDataIOTest
  )

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/AsciiDataIOTest-CorrectWeightedEmom.dat 
    ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiDataIOTest-CorrectWeightedEmom.dat 
    COPYONLY
  )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/AsciiDataIOTest-CorrectWeightedMomE.dat 
    ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiDataIOTest-CorrectWeightedMomE.dat 
    COPYONLY
  )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/AsciiDataIOTest-CorrectUnweighted.dat 
    ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiDataIOTest-CorrectUnweighted.dat 
    COPYONLY
  )

endif()