if(TARGET RootData AND TARGET HelicityFormalism)

  set(testSrc FitTest.cpp)

  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  set(testName "Examples_${fileName}")

  #Add compile target
  add_executable(${testName} ${testSrc})

  # Link to Boost libraries AND your targets and dependencies
  target_link_libraries(
    ${testName}
    PRIVATE Core
            Minuit2IF
            MinLogLH
            RootData
            Integration
            HelicityFormalism
            Boost::unit_test_framework)

  target_include_directories(${testName} PUBLIC ${Boost_INCLUDE_DIR})

  # Move testing binaries into a testBin directory
  set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                               ${PROJECT_BINARY_DIR}/bin/test/)

  # Finally add it to test execution -
  # Notice the WORKING_DIRECTORY and COMMAND
  add_test(
    NAME ${testName}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName})

  configure_file(${PROJECT_SOURCE_DIR}/Physics/particle_list.xml
                 ${PROJECT_BINARY_DIR}/bin/test/particle_list.xml COPYONLY)

else()
  message(WARNING "Requirements not found! Not building Examples Test module!")
endif()
