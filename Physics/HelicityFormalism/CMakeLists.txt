##############################
# Helicity Formalism library #
##############################

#
# Find Boost
#
find_package( Boost COMPONENTS
               unit_test_framework
               system
               REQUIRED QUIET
)

SET( lib_srcs IncoherentIntensity.cpp CoherentIntensity.cpp
    AmpWignerD.cpp PartialDecay.cpp SequentialTwoBodyDecay.cpp
    HelicityKinematics.cpp )

SET( lib_headers IncoherentIntensity.hpp CoherentIntensity.hpp
    SequentialTwoBodyDecay.hpp PartialDecay.hpp PartialDecay.hpp
    SequentialTwoBodyDecay.hpp AmpWignerD.hpp HelicityKinematics.hpp
    ../Amplitude.hpp ../Resonance.hpp )

add_library( HelicityFormalism
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries( HelicityFormalism
  Core
  qft++
  DecayDynamics
  ${Boost_LIBRARIES}
)

target_include_directories( HelicityFormalism
  PUBLIC ${Boost_INCLUDE_DIR}
)

install (FILES ${lib_headers}
    DESTINATION include/ComPWA/Physics/HelicityFormalism
)

install(TARGETS HelicityFormalism
        LIBRARY DESTINATION lib
)

# ------ TESTING ----- #
#Testing routines are stored in separate directory
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)


# Check if Tools target is going to be build. Skip testing section otherwise.
IF( TARGET Tools )

  #Run through each source
  foreach(testSrc ${TEST_SRCS})

    #Extract the filename without an extension (NAME_WE)
    get_filename_component(fileName ${testSrc} NAME_WE)
    SET(testName "HelicityFormalismTest_${fileName}")

    #Add compile target
    add_executable( ${testName} ${testSrc} )


    #link to Boost libraries AND your targets and dependencies
    target_link_libraries( ${testName}
      Core
      HelicityFormalism
      Tools
      pthread
      ${Boost_LIBRARIES}
      ${ROOT_LIBRARIES}
    )
  
    target_include_directories( ${testName}
      PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

    #I like to move testing binaries into a testBin directory
    set_target_properties( ${testName}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/testBin
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
        ${EXECUTABLE_OUTPUT_PATH}/testBin/${TestInName} )
    endforeach( TestIn )

    #Finally add it to test execution -
    #Notice the WORKING_DIRECTORY and COMMAND
    add_test(NAME ${testName}
      WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/testBin
      COMMAND ${EXECUTABLE_OUTPUT_PATH}/testBin/${testName} )

  endforeach(testSrc)

ELSE ()
  MESSAGE( WARNING "Required targets not found! Not building Helicity\
                   formalism testing module!")
ENDIF()
