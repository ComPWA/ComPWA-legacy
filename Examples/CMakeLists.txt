add_subdirectory( DalitzAnalysis_D0KsKK )
add_subdirectory( DalitzFit )
add_subdirectory( Benchmark )
add_subdirectory( SimFit )
add_subdirectory( CanonicalSum )
add_subdirectory( EvtGen )
#add_subdirectory( SliceFit )

if( ${ROOT_FOUND} AND ${GSL_FOUND}
    AND TARGET DataReader AND TARGET HelicityFormalism
)

# ------ TESTING ----- #
#Testing routines are stored in separate directory
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

if(${CMAKE_HOST_APPLE})
	link_directories(${ROOT_LIBRARY_DIR})
endif()

#Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  set(testName "Examples_${fileName}")

  #Add compile target
  add_executable( ${testName} ${testSrc} )

  # Link to Boost libraries AND your targets and dependencies
  target_link_libraries( ${testName}
    Core
    Minuit2IF
    MinLogLH
    Tools
    HelicityFormalism
    Boost::unit_test_framework
  )
  
  target_include_directories( ${testName}
    PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR} )

  # Move testing binaries into a testBin directory
  set_target_properties( ${testName}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )
  target_compile_definitions(${testName} PRIVATE "BOOST_TEST_DYN_LINK=1")
  
  # Copy input files for test programs - we assume they have the name
  # ${testName}-input*. Multiple files can be copied.
  file(GLOB TestInput
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/${fileName}-input*
  )
  foreach( TestIn ${TestInput} )
    get_filename_component( TestInName ${TestIn} NAME )

    add_custom_command(
      TARGET ${testName} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/test/${TestInName}
      ${PROJECT_BINARY_DIR}/bin/test/${TestInName} )
  endforeach( TestIn )

  # Finally add it to test execution -
  # Notice the WORKING_DIRECTORY and COMMAND
  add_test(NAME ${testName}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName} )

  # Do not install testing binaries
  #INSTALL(TARGETS ${testName}
  #  RUNTIME DESTINATION bin/test
  #)
endforeach(testSrc)

else ()
  message(WARNING "Requirements not found! Not building Examples Test module!")
endif()
