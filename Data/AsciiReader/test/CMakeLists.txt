if(TARGET AsciiReader)

  add_executable(Data_AsciiReaderTest AsciiReaderTest.cpp)

  target_link_libraries(Data_AsciiReaderTest
    AsciiReader EvtGenGenerator
    Boost::unit_test_framework
  )
  set_target_properties(Data_AsciiReaderTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(NAME Data_AsciiReaderTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiReaderTest
  )

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/AsciiReaderTest-CorrectWeightedEmom.dat 
    ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiReaderTest-CorrectWeightedEmom.dat 
    COPYONLY
  )
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/AsciiReaderTest-CorrectWeightedMomE.dat 
    ${PROJECT_BINARY_DIR}/bin/test/Data_AsciiReaderTest-CorrectWeightedMomE.dat 
    COPYONLY
  )

endif()