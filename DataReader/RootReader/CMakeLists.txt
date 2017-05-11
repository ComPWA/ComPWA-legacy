# Create RootReader library.
find_package( ROOT REQUIRED )
MESSAGE( STATUS "FindROOT libraries: " "${ROOT_LIBRARIES}/" )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_library (ComPWA_RootReader SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_RootReader ComPWA_DataReader ${ROOT_LIBRARIES} )

target_include_directories (ComPWA_RootReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR})

#install (TARGETS ComPWA_RootReader DESTINATION lib/ComPWA)
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
  target_link_libraries( ${testName} ComPWA_Core ComPWA_RootReader ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )

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
