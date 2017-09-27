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
# Install libraries
#
INSTALL(FILES ${lib_headers} DESTINATION include/ComPWA)

INSTALL(TARGETS Core
        LIBRARY DESTINATION lib
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
  #INSTALL(TARGETS ${testName}
  #  RUNTIME DESTINATION bin
  #)

  #link to Boost libraries AND your targets and dependencies
  TARGET_LINK_LIBRARIES( ${testName}
    Core
    pthread
    ${Boost_LIBRARIES}
  )

  #I like to move testing binaries into a testBin directory
  set_target_properties( ${testName}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/testBin
  )
  set_target_properties( ${testName}
    PROPERTIES PREFIX ""
  )

  # Copy input files for test programs - we assume they have the name
  # ${testName}-input*. Multiple files can be copied.
  FILE(GLOB TestInput
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

#FIND_PACKAGE( SWIG )
#IF(SWIG_FOUND)
#    INCLUDE(${SWIG_USE_FILE})
#
#    FIND_PACKAGE(PythonLibs)
#    IF(PYTHONLIBS_FOUND)
#        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
#
#        INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#
#        SET(CMAKE_SWIG_FLAGS "")
#
#        SET_SOURCE_FILES_PROPERTIES(CoreToScript.i PROPERTIES CPLUSPLUS ON)
#        #SET_SOURCE_FILES_PROPERTIES(CoreToScript.i PROPERTIES SWIG_FLAGS "-includeall")
#        SWIG_ADD_MODULE(core python CoreToScript.i)
#        SWIG_LINK_LIBRARIES(core ${PYTHON_LIBRARIES})
#    ENDIF()
#ENDIF()

#IF(PYTHONLIBS_FOUND)
    # Define the wrapper library that wraps our library
#    ADD_LIBRARY( Core_ext SHARED Core_ext.cpp )
#    TARGET_LINK_LIBRARIES( Core_ext ${Boost_LIBRARIES} Core )
    # don't prepend wrapper library name with lib
#    set_target_properties( Core_ext PROPERTIES PREFIX "" )
#ENDIF()
