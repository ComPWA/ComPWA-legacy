################################
# Create DecayDynamics library #
################################

SET( lib_srcs AbstractDynamicalFunction.cpp RelativisticBreitWigner.cpp
    NonResonant.cpp AmpFlatteRes.cpp Voigtian.cpp Utils/Faddeeva.cc)

SET( lib_headers AbstractDynamicalFunction.hpp
    NonResonant.hpp RelativisticBreitWigner.hpp
    AmpFlatteRes.hpp Voigtian.hpp Utils/Faddeeva.hh FormFactor.hpp )

add_library( DecayDynamics
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries( DecayDynamics
  Core
)

target_include_directories( DecayDynamics
  PUBLIC ${Boost_INCLUDE_DIR}
)

#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/ComPWA/Physics/DecayDynamics
)

install(TARGETS DecayDynamics
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
  SET(testName "HelicityFormalismTest_${fileName}")

  # Add compile target
  add_executable( ${testName} ${testSrc})

  #link to Boost libraries AND your targets and dependencies
  target_link_libraries( ${testName}
    Core
    HelicityFormalism
    pthread
    ${Boost_LIBRARIES}
  )

  # Move testing binaries into a testBin directory
  set_target_properties( ${testName}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )

  # Copy input files for test programs - we assume they have the name
  # ${fileName}-input*. Multiple files can be copied.
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

  #Finally add it to test execution -
  #Notice the WORKING_DIRECTORY and COMMAND
  add_test(NAME ${testName}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName} )
endforeach(testSrc)
