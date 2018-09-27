if(TARGET Tools)

if(${CMAKE_HOST_APPLE})
	link_directories(${ROOT_LIBRARY_DIR})
endif()

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
  # Run through each source
  foreach(testSrc ${TEST_SRCS})

    # Extract the filename without an extension (NAME_WE)
    get_filename_component(fileName ${testSrc} NAME_WE)
    set(testName "HelicityFormalismTest_${fileName}")

    # Add compile target
    add_executable(${testName} ${testSrc})

    # Link to Boost libraries AND your targets and dependencies
    target_link_libraries(${testName}
      Core
      HelicityFormalism
      Tools
      Boost::unit_test_framework
    )
  
    target_include_directories(${testName}
      PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

    # Move testing binaries into a testBin directory
    set_target_properties( ${testName}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    )
    target_compile_definitions(${testName} PRIVATE "BOOST_TEST_DYN_LINK=1")

    # Copy input files for test programs - we assume they have the name
    # ${fileName}-input*. Multiple files can be copied.
    file(GLOB TestInput
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${fileName}-input*
    )
    foreach( TestIn ${TestInput} )
      get_filename_component( TestInName ${TestIn} NAME )

      add_custom_command(
        TARGET ${testName} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/${TestInName}
        ${PROJECT_BINARY_DIR}/bin/test/${TestInName} )
    endforeach( TestIn )

    # Finally add it to test execution -
    # Notice the WORKING_DIRECTORY and COMMAND
    add_test(NAME ${testName}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
      COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName} )

  endforeach(testSrc)

else()
  message(WARNING "Requirements not found! Not building tests!")
endif()