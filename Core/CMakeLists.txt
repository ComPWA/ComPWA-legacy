#######################
# ComPWA Core library #
#######################

#
# FindBoost. Add 'log' dependency. Note
# that all dependencies have to be listed,
# even if they already appear in global
# CMakeLists.txt.
#
FIND_PACKAGE( Boost COMPONENTS
               unit_test_framework
               serialization
               log
               system
               date_time
               REQUIRED QUIET
)

#
# Create library
#
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

ADD_LIBRARY( Core
  SHARED ${lib_srcs} ${lib_headers}
)

TARGET_LINK_LIBRARIES( Core
  ${Boost_LIBRARIES}
)

TARGET_INCLUDE_DIRECTORIES( Core
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR}
)

#
# Install
#
INSTALL(FILES ${lib_headers}
    DESTINATION include/ComPWA/Core
)

INSTALL(TARGETS Core
        LIBRARY DESTINATION lib/ComPWA
)

#
# TESTING
#
#Testing routines are stored in separate directory
FILE(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  SET(testName "CoreTest_${fileName}")

  #Add compile target
  add_executable( ${testName} ${testSrc})

  # Link to Boost libraries AND your targets and dependencies
  TARGET_LINK_LIBRARIES( ${testName}
    Core
    pthread
    ${Boost_LIBRARIES}
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
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName} )

  # Do not install testing binaries
  #INSTALL(TARGETS ${testName}
  #  RUNTIME DESTINATION bin/test
  #)
endforeach(testSrc)
