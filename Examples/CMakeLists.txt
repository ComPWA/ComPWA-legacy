add_subdirectory( DalitzAnalysis_D0KsKK )
add_subdirectory( DalitzFit )
add_subdirectory( Benchmark )
add_subdirectory( SimFit )
add_subdirectory( CanonicalSum )
add_subdirectory( EvtGen )
#add_subdirectory( SliceFit )

#if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#  add_subdirectory( PythonDalitzFit )
#endif()

#
# Find Boost
#
find_package( Boost COMPONENTS
               unit_test_framework
               serialization
               REQUIRED QUIET
)

#
# Find ROOT
#
find_package( ROOT QUIET )

IF( ${ROOT_FOUND} AND ${GSL_FOUND}
    AND TARGET DataReader AND TARGET HelicityFormalism
)

# ------ TESTING ----- #
#Testing routines are stored in separate directory
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  SET(testName "Examples_${fileName}")

  #Add compile target
  add_executable( ${testName} ${testSrc} )

  # Link to Boost libraries AND your targets and dependencies
  TARGET_LINK_LIBRARIES( ${testName}
    Core
    Minuit2IF
    MinLogLH
    Tools
    HelicityFormalism
    ${ROOT_LIBRARIES}
	-pthread
    ${Boost_LIBRARIES}
  )
  
  target_include_directories( ${testName}
    PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR} )

  # Move testing binaries into a testBin directory
  set_target_properties( ${testName}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )
  
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

ELSE ()
  MESSAGE(WARNING "Requirements not found! Not building Examples Test module!")
ENDIF()
