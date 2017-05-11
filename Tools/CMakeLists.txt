# Create DalitzKinematics library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)


### Generator dictionaries
# Currently we do not generate dictionaries of compatibility reasons
# We have to add this directory to the include directories
# so rootcint can find the header and LinkDef files
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#ROOT_GENERATE_DICTIONARY( G__Dict SimpleEfficiency.hpp LINKDEF libDPKinematicsLinkDef.h )
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/G__Dict.cxx.o DESTINATION lib/)
#add_library ( DPKinematics SHARED ${lib_srcs} G__Dict.cxx ${lib_headers} )

add_library ( Tools
  SHARED ${lib_srcs} ${lib_headers} )

target_link_libraries( Tools
  Core
  DataReader
  RootReader
  ${ROOT_LIBRARIES}
  ${Boost_LIBRARIES}
  ${GSL_LIBRARIES}
)

target_include_directories( Tools
  PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR} )


install (FILES ${lib_headers} DESTINATION include/ComPWA)
install(TARGETS Tools
        LIBRARY DESTINATION lib
)


# ------ TESTING ----- #
#Testing routines are stored in separate directory
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

#Run through each source
foreach(testSrc ${TEST_SRCS})
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(fileName ${testSrc} NAME_WE)
  SET(testName "${fileName}")

  #Add compile target
  add_executable( ${testName} ${testSrc})
  #install(TARGETS ${testName}
  #  RUNTIME DESTINATION bin
  #)

  #link to Boost libraries AND your targets and dependencies
  target_link_libraries( ${testName}
    Core
    Tools
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
