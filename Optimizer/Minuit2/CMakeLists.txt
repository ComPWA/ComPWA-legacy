####################################
# Create Minuit2 interface library #
####################################

#
# FindMinuit2
#
find_package( Minuit2 REQUIRED QUIET )

#
# Find Boost
#
find_package( Boost COMPONENTS
               unit_test_framework
               system
               serialization
               REQUIRED QUIET
)

set( lib_srcs
     MinuitFcn.cpp MinuitIF.cpp MinuitResult.cpp
)

set( lib_headers
     MinuitFcn.hpp MinuitIF.hpp MinuitResult.hpp ../Optimizer.hpp
)

add_library( Minuit2IF
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries( Minuit2IF
  Core
  ${Boost_LIBRARIES}
  ${MINUIT2_LIBRARIES}
)

target_include_directories( Minuit2IF
  PUBLIC ${MINUIT2_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)

#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/ComPWA/Optimizer/Minuit2IF
)

install (TARGETS Minuit2IF
    LIBRARY DESTINATION lib/ComPWA
)

#
# TESTING
#
# Testing routines are stored in separate directory
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

# Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  SET(testName "OptimizerTest_${fileName}")

  # Add compile target
  add_executable( ${testName} ${testSrc})

  # Link to Boost libraries AND your targets and dependencies
  target_link_libraries( ${testName}
    Core
    RootReader
    ${Boost_LIBRARIES}
    ${MINUIT2_LIBRARIES}
  )

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
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName} )
endforeach(testSrc)
