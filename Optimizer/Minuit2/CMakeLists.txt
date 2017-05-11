# Create Minuit2IF library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED COMPONENTS Minuit2 )
find_package( GSL REQUIRED gsl gslcblas )

SET(lib_srcs MinuitFcn.cpp MinuitIF.cpp MinuitResult.cpp)
SET(lib_headers MinuitFcn.hpp MinuitIF.hpp MinuitResult.hpp ../Optimizer.hpp)

add_library( Minuit2IF
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries( Minuit2IF
  Core
  ${Boost_LIBRARIES}
  ${ROOT_LIBRARIES}
  ${GSL_LIBRARIES}
)

target_include_directories( Minuit2IF
  PUBLIC ${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)

# ------ TESTING ----- #
#Testing routines are stored in separate directory
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  SET(testName "OptimizerTest_${fileName}")

  #Add compile target
  add_executable( ${testName} ${testSrc})

  #link to Boost libraries AND your targets and dependencies
  target_link_libraries( ${testName}
    Core
    RootReader
    ${Boost_LIBRARIES}
    ${ROOT_LIBRARIES}
  )

  #I like to move testing binaries into a testBin directory
  set_target_properties( ${testName}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/testBin
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
      ${EXECUTABLE_OUTPUT_PATH}/testBin/${TestInName} )
  endforeach( TestIn )

  #Finally add it to test execution -
  #Notice the WORKING_DIRECTORY and COMMAND
  add_test(NAME ${testName}
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/testBin
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/testBin/${testName} )
endforeach(testSrc)
