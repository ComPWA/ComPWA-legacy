if(TARGET RootData AND TARGET HelicityFormalism)

  set(TEST_SRC FitTest.cpp)

  #Extract the filename without an extension (NAME_WE)
  get_filename_component(
    FILENAME
    ${TEST_SRC}
    NAME_WE
  )
  set(TEST_NAME "Examples_${FILENAME}")

  # Add compile target
  add_executable(${TEST_NAME} ${TEST_SRC})

  # Link to Boost libraries AND your targets and dependencies
  target_link_libraries(
    ${TEST_NAME}
    PRIVATE Core
            Minuit2IF
            MinLogLH
            RootData
            Integration
            HelicityFormalism
            Boost::unit_test_framework
  )

  target_include_directories(${TEST_NAME} PUBLIC ${Boost_INCLUDE_DIR})

  # Move testing binaries into a testBin directory
  set_target_properties(
    ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                            ${PROJECT_BINARY_DIR}/bin/test/
  )

  # Finally add it to test execution -
  # Notice the WORKING_DIRECTORY and COMMAND
  add_test(
    NAME ${TEST_NAME}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${TEST_NAME}
  )

  configure_file(
    ${PROJECT_SOURCE_DIR}/Physics/particle_list.xml
    ${PROJECT_BINARY_DIR}/bin/test/particle_list.xml
    COPYONLY
  )

else()
  message(WARNING "Requirements not found! Not building Examples Test module!")
endif()
