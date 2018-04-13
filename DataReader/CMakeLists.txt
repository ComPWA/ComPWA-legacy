#############################
# Create DataReader library #
#############################

SET(lib_srcs
  Data.cpp
  DataCorrection.cpp
  CorrectionTable.cpp
)
SET(lib_headers
  Data.hpp
  DataCorrection.hpp
  CorrectionTable.hpp
)

add_library( DataReader
  SHARED ${lib_srcs} ${lib_headers}
)

#
# AsciiReader module
#
add_subdirectory (AsciiReader)

#
# RootReader module
#
find_package( ROOT QUIET )
IF( ${ROOT_FOUND} )
    add_subdirectory (RootReader)
ELSE ()
  MESSAGE(WARNING "ROOT not found! Not building RootReader module!")
ENDIF()


target_link_libraries( DataReader
  Core ${Boost_LIBRARIES}
)

target_include_directories( DataReader
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/ComPWA/DataReader
)

install(TARGETS DataReader
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
  SET(testName "${fileName}")

  # Add compile target
  add_executable( ${testName} ${testSrc})

  #link to Boost libraries AND your targets and dependencies
  target_link_libraries( ${testName}
    Core
    RootReader
    ${Boost_LIBRARIES}
  )

  # Move testing binaries into a testBin directory
  set_target_properties( ${testName}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
  )

  # Copy input files for test programs - we assume they have the name
  # ${testName}-input*. Multiple files can be copied.
  file(GLOB TestInput RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/${testName}-input*)
  foreach( TestIn ${TestInput} )
    get_filename_component( TestInName ${TestIn} NAME )

    add_custom_command(
      TARGET ${testName} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/test/${TestInName}
      ${PROJECT_BINARY_DIR}/bin/test/${TestInName} )
  endforeach( TestIn )

  #Finally add it to test execution -
  #Notice the WORKING_DIRECTORY and COMMAND
  add_test(NAME ${testName}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName} )
endforeach(testSrc)
