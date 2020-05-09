file(
  GLOB
  SOURCE_LIST
  *.cpp
)
file(
  GLOB
  HEADER_LIST
  *.hpp
)

add_library(FunctionTree ${SOURCE_LIST} ${HEADER_LIST})

target_include_directories(FunctionTree PUBLIC Boost::serialization)

target_link_libraries(
  FunctionTree
  PUBLIC Core
         pstl::ParallelSTL
         Boost::serialization
)

install(
  TARGETS FunctionTree
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/Core/FunctionTree
)

#
# TESTING
#
# Testing routines are stored in separate directory
file(
  GLOB TEST_SRCS
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  test/*.cpp
)

# Run through each source
foreach(TEST_SRC ${TEST_SRCS})
  # Extract the filename without an extension (NAME_WE)
  get_filename_component(
    FILENAME
    ${TEST_SRC}
    NAME_WE
  )
  set(TEST_NAME "CoreTest_${FILENAME}")

  # Add compile target
  add_executable(${TEST_NAME} ${TEST_SRC})

  # Link to Boost libraries AND your targets and dependencies
  target_link_libraries(
    ${TEST_NAME}
    Core
    FunctionTree
    Threads::Threads
    Boost::unit_test_framework
  )

  target_include_directories(${TEST_NAME} PRIVATE Boost::unit_test_framework)

  # Move testing binaries into a testBin directory
  set_target_properties(
    ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                            ${PROJECT_BINARY_DIR}/bin/test/
  )

  # Copy input files for test programs - we assume they have the name
  # ${TEST_NAME}-input*. Multiple files can be copied.
  file(
    GLOB TEST_INPUT
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    test/${FILENAME}-input*
  )
  foreach(TEST_IN ${TEST_INPUT})
    get_filename_component(
      TEST_IN_NAME
      ${TEST_IN}
      NAME
    )

    add_custom_command(
      TARGET ${TEST_NAME}
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/test/${TEST_IN_NAME}
        ${PROJECT_BINARY_DIR}/bin/test/${TEST_IN_NAME}
    )
  endforeach(TEST_IN)

  # Finally add it to test execution -
  # Notice the WORKING_DIRECTORY and COMMAND
  add_test(
    NAME ${TEST_NAME}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${TEST_NAME}
  )

endforeach(TEST_SRC)
