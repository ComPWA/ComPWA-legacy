# Check if all requirements are found
if(TARGET RootData)

  add_executable(Data_RootDataIOTest RootDataIOTest.cpp)

  target_link_libraries(
    Data_RootDataIOTest
    Boost::unit_test_framework
    ROOT::EG
    ROOT::Physics
    ROOT::RIO
    ROOT::Tree
    RootData
  )

  set_target_properties(
    Data_RootDataIOTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                   ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(
    NAME Data_RootDataIOTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/Data_RootDataIOTest
  )

  add_executable(RootGeneratorTest RootGeneratorTest.cpp)
  target_link_libraries(
    RootGeneratorTest
    Boost::unit_test_framework
    RootData
  )
  set_target_properties(
    RootGeneratorTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                 ${PROJECT_BINARY_DIR}/bin/test/
  )
  add_test(
    NAME RootGeneratorTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/RootGeneratorTest
  )
endif()
