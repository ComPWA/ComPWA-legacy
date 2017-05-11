# Create ComPWA Core library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_library ( ComPWA_Core SHARED ${lib_srcs} ${lib_headers} )
target_link_libraries( ComPWA_Core ${Boost_LIBRARIES} )
target_include_directories ( ComPWA_Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install (TARGETS ComPWA_Core DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)


# ------ TESTING ----- #
#Testing routines are stored in separate directory
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(testName ${testSrc} NAME_WE)

  #Add compile target
  add_executable(${testName} ${testSrc})

  #link to Boost libraries AND your targets and dependencies
  target_link_libraries( ${testName} ${Boost_LIBRARIES} ComPWA_Core pthread )

  #I like to move testing binaries into a testBin directory
  set_target_properties(${testName} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY  ${EXECUTABLE_OUTPUT_PATH}/testBin)

  # Copy input files for test programs - we assume they have the name
  # ${testName}-input*. Multiple files can be copied.
  file(GLOB TestInput RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/${testName}-input*)
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

#find_package( SWIG )
#IF(SWIG_FOUND)
#    INCLUDE(${SWIG_USE_FILE})
#
#    find_package(PythonLibs)
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
#    add_library( Core_ext SHARED Core_ext.cpp )
#    target_link_libraries( Core_ext ${Boost_LIBRARIES} ComPWA_Core )
    # don't prepend wrapper library name with lib
#    set_target_properties( Core_ext PROPERTIES PREFIX "" )
#ENDIF()
