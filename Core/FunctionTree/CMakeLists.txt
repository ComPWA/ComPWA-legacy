file(GLOB lib_srcs *.cpp)
file(GLOB lib_headers *.hpp)

add_library(FunctionTree ${lib_srcs} ${lib_headers})

target_include_directories(FunctionTree PUBLIC Boost::serialization)

target_link_libraries(FunctionTree PUBLIC Core pstl::ParallelSTL
                                          Boost::serialization)

install(
  TARGETS FunctionTree
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/Core/FunctionTree)

#
# TESTING
#
# Testing routines are stored in separate directory
file(
  GLOB TEST_SRCS
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  test/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  set(testName "CoreTest_${fileName}")

  #Add compile target
  add_executable(${testName} ${testSrc})

  # Link to Boost libraries AND your targets and dependencies
  target_link_libraries(${testName} Core FunctionTree Threads::Threads
                        Boost::unit_test_framework)

  target_include_directories(${testName} PRIVATE Boost::unit_test_framework)

  # Move testing binaries into a testBin directory
  set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                               ${PROJECT_BINARY_DIR}/bin/test/)

  # Copy input files for test programs - we assume they have the name
  # ${testName}-input*. Multiple files can be copied.
  file(
    GLOB TestInput
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    test/${fileName}-input*)
  foreach(TestIn ${TestInput})
    get_filename_component(TestInName ${TestIn} NAME)

    add_custom_command(
      TARGET ${testName}
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/test/${TestInName}
        ${PROJECT_BINARY_DIR}/bin/test/${TestInName})
  endforeach(TestIn)

  # Finally add it to test execution -
  # Notice the WORKING_DIRECTORY and COMMAND
  add_test(
    NAME ${testName}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName})

endforeach(testSrc)
