# Check if all requirements are found
if(TARGET RootData)

  add_executable(Data_RootDataIOTest RootDataIOTest.cpp)

  target_link_libraries(
    Data_RootDataIOTest
    RootData
    Boost::unit_test_framework
    ROOT::Physics
    ROOT::EG
    ROOT::Tree
    ROOT::RIO
  )

  set_target_properties(
    Data_RootDataIOTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY
               ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  )

  add_test(
    NAME Data_RootDataIOTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
    COMMAND
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/Data_RootDataIOTest
  )

  add_executable(RootGeneratorTest RootGeneratorTest.cpp)
  target_link_libraries(
    RootGeneratorTest
    RootData
    Boost::unit_test_framework
  )
  set_target_properties(
    RootGeneratorTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY
               ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  )
  add_test(
    NAME RootGeneratorTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
    COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/RootGeneratorTest
  )
endif()
