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
  target_link_libraries( ${testName} ${Boost_LIBRARIES} ComPWA_Core )

  #I like to move testing binaries into a testBin directory
  set_target_properties(${testName} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY  ${EXECUTABLE_OUTPUT_PATH}/testBin)

  #Finally add it to test execution -
  #Notice the WORKING_DIRECTORY and COMMAND
  add_test(NAME ${testName}
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/testBin
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/testBin/${testName} )
endforeach(testSrc)
